volcan = function(res, min.lfc=0, max.pval=0.05, lim.pval=0, lim.lfc=c(-Inf, Inf), label.names=c(''), col.insign='#99999999', col.up='#cc5544', col.down='dodgerblue3', main='', use.adj.pval=T) {
   if(use.adj.pval) { res$p_val = res$p_val_adj }
   if(length(min.lfc)==1) min.lfc=c(-min.lfc, min.lfc)
   if(lim.lfc[1]==-Inf) lim.lfc[1] = min(res$avg_log2FC)
   if(lim.lfc[2]==Inf) lim.lfc[2] = max(res$avg_log2FC)
   res$sym = 'unchanged'
   res$Gene = row.names(res)
   res = within(res, {
      reg = factor(rep('insign', nrow(res)), levels=c('down', 'insign', 'up'))
      reg[p_val <= max.pval & avg_log2FC <= min.lfc[1]] = 'down'
      reg[p_val <= max.pval & avg_log2FC >= min.lfc[2]] = 'up'
      if(is.logical(label.names)) {
         if(label.names) {
            Gene[reg=='insign'] = NA
         } else {
            Gene = NA
         }
      } else {
         Gene[!Gene %in% label.names] = NA
      }
      sym[avg_log2FC < lim.lfc[1]] = 'left'
      sym[avg_log2FC > lim.lfc[2]] = 'right'
      sym[p_val < lim.pval] = 'top'
      sym[p_val < lim.pval & avg_log2FC < lim.lfc[1]] = 'topleft'
      sym[p_val < lim.pval & avg_log2FC > lim.lfc[2]] = 'topright'
      sym = factor(sym)
      avg_log2FC[avg_log2FC < lim.lfc[1]] = lim.lfc[1]
      avg_log2FC[avg_log2FC > lim.lfc[2]] = lim.lfc[2]
      p_val[p_val < lim.pval] = lim.pval
   })
   lab = data.frame(x=lim.lfc, y=c(1, 1), label=as.character(c(sum(res$reg=='down'), sum(res$reg=='up'))), reg=c('down', 'up'))
   return(ggplot(data=res, aes(x=avg_log2FC, y=-log10(p_val), col=reg, label=Gene, shape=sym)) + geom_point(size=2) + theme_classic() +
             geom_text_repel(data=res, aes(x=avg_log2FC, y=-log10(p_val), label=Gene), point.size=3, show.legend=F, na.rm=T, size=3, point.padding=0.2, box.padding=1, segment.curvature=-1e-20, max.overlaps=Inf,
                             force=1.5, force_pull=0.4, segment.size=0.1, min.segment.length=unit(0, 'lines'), nudge_x=ifelse(res$reg=='down', -1, 1), nudge_y=ifelse(res$reg=='insign', -0.1, .1)) +
             annotate(geom='text', x=lim.lfc, y=c(0,0), label=c(sum(res$reg=='down'), sum(res$reg=='up')), col=c(col.down, col.up)) +
             scale_color_manual(values=list(up=col.up, down=col.down, insign=col.insign)) + scale_y_continuous(expand=c(0,0.1)) +
             scale_shape_manual(values=c("unchanged"="\u25CF", "top"="\u25B2", "left"="\u25C0", "right"="\u25B6", "topleft"="\u25E4", "topright"="\u25E5"), guide='none') +
             ggtitle(main, subtitle=ifelse(use.adj.pval, '* adjusted p-values', '* un-adjusted p-values')) + geom_vline(xintercept=min.lfc, lty=2, color=col.insign) + geom_hline(yintercept=-log10(max.pval), lty=2, color=col.insign))
}


vennDE = function(results, max_p=0.05, min_log2FC=0, use.adj.pval=F, title='overlap of differentially expressed genes') {
   for(i in names(results)) {
      results[[i]] = within(results[[i]], {p_val[is.na(p_val)] = 1; p_val_adj[is.na(p_val_adj)] = 1; if(use.adj.pval) p_val = p_val_adj})
   }
   # filter reg, up, down from results per comparison
   reg = lapply(results, function(r) {list(reg=row.names(r[r$p_val<max_p,]), up=row.names(r[r$p_val<max_p & r$avg_log2FC>min_log2FC,]), down=row.names(r[r$p_val<max_p & r$avg_log2FC<(-min_log2FC),]))})
   print(ggvenn(lapply(reg, function(r) {r$reg}), set_name_size=4, stroke_size=0.5, show_percentage=F, fill_color=as.vector(cols$Health)) + ggtitle(title))
   
   comps = c()
   for(i in names(reg)) {
      other = names(reg)[names(reg) != i] # other dataset names
      cat(paste(i, 'exclusive:   up-regulated:', paste(reg[[i]]$up[!reg[[i]]$up %in% unlist(lapply(other, function(o) reg[[o]]$reg))], collapse=' '),
                '  down-regulated:', paste(reg[[i]]$down[!reg[[i]]$down %in% unlist(lapply(other, function(o) reg[[o]]$reg))], collapse=' '), '\n'))
      wrong = reg[[i]]$up[ reg[[i]]$up %in% reg[[other[1]]]$down & reg[1]$up %in% reg[other[2]]$down ]
      for (w in wrong) warning(paste('there is conflicting regulation for', w, ': up in', i, 'down in the others'))
      wrong = reg[[i]]$down[ reg[[i]]$down %in% reg[[other[1]]]$up & reg[[i]]$down %in% reg[[other[2]]]$up]
      for (w in wrong) warning(paste('there is conflicting regulation for', w, ': down in', i, 'up in the others'))
      for(o in other) {
         c = paste(sort(c(i, o)), collapse=' & ')
         if(c %in% comps) next # already compared
         wrong = c(reg[[i]]$up[reg[[i]]$up %in% reg[[o]]$down], reg[[i]]$down[reg[[i]]$down %in% reg[[o]]$up])
         for (w in wrong) warning(paste('there is conflicting regulation for', w, 'between', c))
         rest = other[other != o]
         if(length(rest)>0) rest = unlist(lapply(rest, function(r) reg[[r]]$reg))
         cat(paste(c, '   up-regulated:', paste(reg[[i]]$up[reg[[i]]$up %in% reg[[o]]$up & !reg[[i]]$up %in% rest], collapse=' '), '   down-regulated:', paste(reg[[i]]$down[reg[[i]]$down %in% reg[[o]]$down & !reg[[i]]$down %in% rest], collapse=' '),'\n'))
         comps = c(comps, c)
      }
   }
   if(length(results)>2) 
      cat(paste('in ALL    up-regulated:', paste(reg[[i]]$up[sapply(reg[[i]]$up, function(x) all(sapply(reg, function(r) x %in% r$up)))], collapse=' '),
                '   down-regulated:', paste(reg[[i]]$down[sapply(reg[[i]]$down, function(x) all(sapply(reg, function(r) x %in% r$down)))], collapse=' '),'\n'))
}


scHeat = function(so.sub, var, first, second, plot=T) {
   so.c.cells = table(so.sub@meta.data$Sample_Name)
   Idents(so.sub) = so.sub@meta.data[[var]]
   if(length(table(Idents(so.sub)))==1 | any(table(Idents(so.sub))[c(first, second)]<3)) { cat(paste(c('\nless than 3 cells per condition\n'))); return() }
   
   so.c.cells = table(so.sub@meta.data$Sample_Name)
   so.c.mark = FindMarkers(so.sub, test.use="MAST", ident.1=first, ident.2=second, assay='RNA', verbose=F, min.pct=0.1)
   so.c.mark = subset(so.c.mark, p_val_adj<0.1 & abs(avg_log2FC)>0.1)
   if(nrow(so.c.mark)<1) { cat(paste(c('\nno significant genes for', first, 'vs', second,'\n'))); return();}
   
   l2fc = so.c.mark[,"avg_log2FC", drop=F]
   padj = so.c.mark[,"p_val_adj", drop=F]
   padj = ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,"*","")))
   
   p1 = Heatmap(l2fc, name="l2FC", column_gap=unit(0, "mm"), border=T, column_names_gp=gpar(fontsize=10),
                row_names_gp=gpar(fontsize=10), column_title_gp=gpar(fontsize=8), row_title_rot=0, 
                cluster_rows=T, cluster_columns=F, show_heatmap_legend=F,
                cell_fun=function(j, i, x, y, w, h, fill) { grid.text(paste0(signif(l2fc[i,j], 2), " ", padj[i,j]), x, y, gp=gpar(fontsize=7)) },
                col=colorRamp2(c(-2, 0, 2), c("#3344aa", "#eeeeee", "#aa4444")))
   
   pct = so.c.mark[,c("pct.1", "pct.2"), drop=F]
   colnames(pct) = c(first, second)
   
   p2 = Heatmap(as.matrix(pct), name="pct", column_names_gp=gpar(fontsize=10),
                row_names_gp=gpar(fontsize=12), column_title_gp=gpar(fontsize=8), column_gap=unit(0, "mm"),
                cluster_rows=T, cluster_columns=F, row_title_rot=0, border=T, show_heatmap_legend=F,
                cell_fun=function(j, i, x, y, w, h, fill) { grid.text(signif(pct[i, j], 2), x, y, gp=gpar(fontsize=10)) },
                col=colorRamp2(c(0, 1), c("#ffffff", "#555555")))
   
   avg_expr = AverageExpression(so.sub, return.seurat=F, assays="RNA", group.by="Sample_Name", slot="data", verbose=F)[[1]][rownames(so.c.mark), , drop=F]
   avg_expr = avg_expr/rowSums(avg_expr, na.rm=T)*100
   
   p3 = Heatmap(as.matrix(avg_expr), name="mean(exp)", row_km=2, row_names_gp=gpar(fontsize=10),
                cluster_rows=T, cluster_columns=F, show_heatmap_legend=F, column_names_gp=gpar(fontsize=10),
                column_title_gp=gpar(fontsize=8), column_gap=unit(0, "mm"), row_title_rot=0, border=T,
                cell_fun=function(j, i, x, y, w, h, fill) { grid.text(paste(round(avg_expr[i, j], 1), '%'), x, y, gp=gpar(fontsize=10)) },
                top_annotation=columnAnnotation(cells=anno_barplot(unlist(list(so.c.cells)), height=unit(10, "mm"), border=F, add_numbers=T, gp=gpar(fill="#CCCCCC", col='white'))),
                col=colorRamp2(c(0, max(avg_expr)), c("#ffffff", "#44aa55")))
   
   if(plot) {print(p1+p2+p3)}
   so.c.mark
}


cols = list(Health=c(PV='#E0ED7C', BP='#2BCCC1', PSORI='#6AC631', `PV&BP`='#70eD9C', `PV&PSORI`='#a0dD5C', `BP&PSORI`='#20CD8C', `PV&BP&PSORI`='#F0FDDC', HC='#C5CAC7'), Lesion=c(lesional='#ff5577', perilesional='#ff99aa', ' '='#555555', healthy='#C5CAC7'), Replicate=c('1'='#ccaaff', '2'='#8866cc', '3'='#440077', '4'='#110055'), dir=c(up='#ff4466', down='#4466ff'))

biramp = function(n) {as.vector(t(cbind(viridis(ceiling(n/2)), colorRampPalette(c('#ddee11', '#ff5533', '#4422bb'), interpolate="spline", space="Lab")(ceiling(n/2)))))}

subchunkify <- function(g, fig_height=7, fig_width=5) {
   if(!isTRUE(getOption('knitr.in.progress'))){ print(g); return(); }
   g_deparsed <- paste0(deparse(function() {g}), collapse='')
   if(!exists("chunky_id")){ chunky_id <<- 0}
   chunky_id <<- chunky_id+1
   sub_chunk <- paste0("\n`","``{r sub_chunky_", chunky_id, ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                       "\n(", g_deparsed, ")()", "\n`","``\n")
   cat(knitr::knit(text = knitr::knit_expand(text=sub_chunk), quiet=T))
}