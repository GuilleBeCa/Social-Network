########## TASK 2.1 Complete Mixing

readjustment = function(x, pair, d, mu){
    i = pair[1]
    j = pair[2]
    xi = x[i]
    xj = x[j]
    if (abs(xi - xj) < d){
        x[i] = xi + mu * (xj - xi)
        x[j] = xj + mu * (xi - xj)
    } else {
        x[i] = xi
        x[j] = xj
    }
    return(x)
}

updating_opinions = function(N, T, npair, d, mu){
    initial_opinions = runif(N, min=0, max=1)
    data = matrix(, nrow=N, ncol=T)
    data[, 1] = initial_opinions
    intermediate_opinions = initial_opinions

    for (t in 2:T){
        for (n in 1:npair){
            pair = sample.int(N, 2, replace=FALSE)
            intermediate_opinions = readjustment(intermediate_opinions, pair, d, mu)
        }
        data[, t] = intermediate_opinions
    }
    df_opinions = data.frame(TimeStep=rep(1:T, each=N), Opinion=c(data))
    return(df_opinions)
}

plot_time_series = function(df_opinions){
    g = ggplot(df_opinions, aes(x=TimeStep, y=Opinion)) +
        geom_point() +
        xlab("Time") +
        scale_x_continuous(limits=c(0, 1.01 * T),
            expand=c(0, 0), 
            breaks=round(seq(0, T, by=T%/%20))) +
        ylab("Opinion") + 
        ylim(c(0, 1)) + 
        ggtitle(paste0("Time series with d = ", d, ", µ = ", mu, ": ", T, " time steps with ", npair, " pair samplings at every time step, for ", N, " agents.")) + 
        theme(axis.text.x=element_text(size=24), axis.text.y=element_text(size=24), plot.title=element_text(size=22), axis.title=element_text(size=24))
    return(g)
}

plot_final_vs_initial_opinions = function(df_opinions){
    T = max(df_opinions$TimeStep)
    initial_opinions = df_opinions[df_opinions$TimeStep == 1,]$Opinion
    final_opinions = df_opinions[df_opinions$TimeStep == T,]$Opinion
    df = data.frame(initial_opinions, final_opinions)
    
    g = ggplot(df, aes(x=initial_opinions, y=final_opinions)) + 
        geom_point() + 
        xlab("Initial opinion") + 
        ylab("Final opinion") + 
        scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.1)) + 
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
        # ggtitle(paste0("Final vs initial opinions with d = ", d, ", µ = ", mu, ": ", npair * (T - 1), " total pair samplings, for ", N, " agents.")) + 
        theme(axis.text.x=element_text(size=34), axis.text.y=element_text(size=34), plot.title=element_text(size=22), axis.title=element_text(size=34))
    return(g)
}

number_of_peaks = function(N, T, npair, d_values, mu, reps, max_peaks){
    colnames = c("d")
    for (i in 1:max_peaks){
        colnames = c(colnames, paste(i, "peak"))
    }
    df_peaks = setNames(data.frame(matrix(ncol=length(colnames), nrow=0)), colnames)

    for (d in d_values){
        set.seed(42)
        peaks_vector = c()
        for (i in 1:reps){
            df_opinions = updating_opinions(N, T, npair, d, mu)
            peaks = df_opinions[df_opinions$TimeStep == T,]$Opinion
            peaks = sort(round(peaks[!duplicated(round(peaks, digits=1))], digits=3))
            peaks_vector = c(peaks_vector, length(peaks))
            # peaks_ch = paste(as.character(peaks), collapse=" ")
            # cat("(", d, ", ", i, "): ", peaks_ch,", ", length(peaks), "\n", sep="")
        }
        v = c(d, table(factor(peaks_vector, levels=1:max_peaks)))
        df_peaks[nrow(df_peaks) + 1,] = v
    }
    if (all(as.numeric(rowSums(df_peaks %>% select(-d))) != reps)){
        return("ERROR!")
    } else {
        return(df_peaks)
    }
}

plot_number_of_peaks_d = function(df_peaks){
    melted = melt(df_peaks, id.vars="d")
    g = ggplot(data=melted, aes(x=d, y=value, color=variable)) + 
        geom_line(linewidth=1) + 
        labs(y="Counts", x="d", color=NULL) +
        scale_x_continuous(limits=c(0, 0.5), breaks=seq(0, 0.5, by=0.05)) + 
        theme(axis.text.x=element_text(size=34), axis.text.y=element_text(size=34), plot.title=element_text(size=22), 
            axis.title=element_text(size=34), legend.title=element_text(size=0), legend.text=element_text(size=22))
    return(g)
}








########## TASK 2.2 Social Networks

sample_pair = function(g){
    edges = E(g)
    random_edge = edges[sample(length(edges), 1)]
    pair = c(ends(g, random_edge))
    return(pair)
}

updating_opinions_lattice = function(g, T, npair, d, mu){
  Nv = vcount(g)
  V(g)$opinions = runif(Nv)
  
  n = as.integer(sqrt(Nv))
  n_list = 1:n
  df_ops = expand.grid(Row=n_list, Col=n_list) # grid of coordinates on lattice
  df_ops["VertexNumber"] = (df_ops$Col - 1) * n + df_ops$Row # increasing number for each vertex 
  df_ops["TimeStep"] = 0
  df_ops["Opinions"] = V(g)$opinions

  pair = sample_pair(g)
  df_ops_t = df_ops[df_ops$TimeStep == 0,]
  updated_opinions = df_ops_t$Opinions

  for (t in 1:T){
    for (i in 1:npair){
      pair = sample_pair(g)
      updated_opinions = readjustment(updated_opinions, pair, d, mu)
    }
    df_ops_t$TimeStep = t
    df_ops_t$Opinions = updated_opinions
    df_ops = rbind(df_ops, df_ops_t)
  }
  return(df_ops)
}

plot_lattice = function(df_ops_t){
    t = unique(df_ops_t$TimeStep)
    p = ggplot(df_ops_t, aes(x=Col, y=Row, fill=Opinions)) +
        geom_tile() + 
        scale_fill_continuous(limits=c(0, 1), breaks=seq(-1, 1, by=0.25)) +
        ggtitle(paste0("Time step: ", t, "/", T)) +
        theme(plot.title=element_text(hjust=0.5, size=38), axis.text.x=element_text(size=34), axis.text.y=element_text(size=34), 
            axis.title=element_text(size=34), legend.title=element_text(size=30), legend.text=element_text(size=26))
    return(p)
}

plot_lattice_grid = function(df_ops, ncol, side, T){
  height = side * (T %/% ncol + 1) 
  width = side * ncol 
  options(repr.plot.width=width, repr.plot.height=height)

  plot_list = list()
  for (t in 0:T){
    df_ops_t = df_ops[df_ops$TimeStep == t,]
    plot_list[[t + 1]] = plot_lattice(df_ops_t)
  }
  title = paste0("Updated opinions on a ", n, " x ", n, " square lattice, with d = ", d, ", µ = ", mu, " and ", npair, " pair samplings per time step")
  return(grid.arrange(grobs=plot_list, ncol=ncol, top=textGrob(title,gp=gpar(fontsize=24))))
}

plot_histogram_final_opinions = function(df_ops){
    T = max(unique(df_ops$TimeStep))
    final_opinions = round(df_ops[df_ops$TimeStep == T,]$Opinions, digits=2)
    max_y = max(as.numeric(table(final_opinions)))
    y_step = max_y%/%10
    if(y_step == 0) {
        y_step = 1
    }

    p = ggplot() + 
        aes(final_opinions) + 
        geom_histogram(binwidth=0.01, fill="blue") +
        xlab("Final opinions") + 
        ylab("Counts") + 
        scale_x_continuous(expand=c(0.05, 0), breaks=seq(0, 1, by=0.05)) + 
        scale_y_continuous(limit=c(0, 1.05 * max_y), expand=c(0, 0), breaks=seq(0, 1.05 * max_y, by=y_step)) + 
        # ggtitle(paste0("Histogram of final opinions for d = ", d, ", µ = ", mu, " on a ", n, " x ", n, " square lattice.")) + 
        theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), plot.title=element_text(size=22), axis.title=element_text(size=20))
    return(p)
}








########## TASK 2.3 Vector opinions

generate_opinion_matrix = function(N, m){
    matrix_opinions = matrix(sample(c(0, 1), N*m, replace=TRUE), nrow=N, ncol=m)
    rownames(matrix_opinions) = 1:N
    colnames(matrix_opinions) = paste0("Subject", 1:m)
    return(matrix_opinions)
}

generate_opinion_df = function(N, m){
    df_opinions = as.data.frame(generate_opinion_matrix(N, m))
    return(df_opinions)
}

# vector_readjustment = function(matrix_opinions, N, m, d, mu){ # one pair at a time
#     pair = sample(1:N, 2, replace=TRUE)
#     selected_pair = matrix_opinions[pair,]
#     agreed_subjects = sum(selected_pair[1,] == selected_pair[2,])

#     if (agreed_subjects > (m - d) & agreed_subjects < m){
#         u = runif(1)
#         if (u < mu){
#             differing_subjects = unname(which(selected_pair[1,] != selected_pair[2,]))
#             if (length(differing_subjects) > 1){
#                 subject = sample(differing_subjects, 1)
#             } else {
#                 subject = differing_subjects
#             }
#             agent = sample(1:2, 1)
#             matrix_opinions[agent, subject] = !matrix_opinions[agent, subject]
#         }
#     }
#     return(matrix_opinions)
# }

vector_readjustment = function(matrix_opinions, N, m, d, mu){ # all pairs of matrix at the same time
    sampling = sample(1:N, N, replace=FALSE)
    reordered_matrix = matrix_opinions[sampling,]
    list = lapply(seq_len(N%/%2),
        function(i){
            pair_readjustment(reordered_matrix[(2*i-1):(2*i),], m, d, mu)
        })

    mo = do.call(rbind, list)
    matrix_opinions = mo[order(as.integer(rownames(mo))),]
    return(matrix_opinions)
}

vector_readjustment_df = function(matrix_opinions, N, m, d, mu){
    df_opinions = as.data.frame(vector_readjustment(matrix_opinions, N, m, d, mu))
    return(df_opinions)
}

pair_readjustment = function(selected_pair, m, d, mu){
    agreed_subjects = sum(selected_pair[1,] == selected_pair[2,])

    if (agreed_subjects >= (m - d) & agreed_subjects < m){
        u = runif(1)
        if (u < mu){
            differing_subjects = unname(which(selected_pair[1,] != selected_pair[2,]))
            if (length(differing_subjects) > 1){
                subject = sample(differing_subjects, 1)
            } else {
                subject = differing_subjects
            }
            agent = sample(1:2, 1)
            selected_pair[agent, subject] = !selected_pair[agent, subject]
        }
    }
    return(selected_pair)
}

generate_final_opinions = function(matrix_opinions, N, m, d, mu, T){
    for (t in 1:T){
        if (nrow(unique(matrix_opinions)) == 1){
            cat("Vector opinions converge at t =", t, "\n")
            break
        }
        matrix_opinions = vector_readjustment(matrix_opinions, N, m, d, mu)
    }
    return(matrix_opinions)
}

generate_final_opinions_df = function(matrix_opinions, N, m, d, mu, T){
    df_opinions = as.data.frame(generate_final_opinions(matrix_opinions))
    return(df_opinions)
}

final_opinions_table = function(matrix_final_opinions){
    final_opinions = apply(unname(matrix_final_opinions), 1, paste, collapse="")
    return(table(final_opinions))
}

distance = function(pair){
    x = as.integer(strsplit(pair[1], split="")[[1]])
    y = as.integer(strsplit(pair[2], split="")[[1]])
    d = sum(x != y)
    return(d)
}

final_opinion_distances = function(final_opinions){
    opinions = names(final_opinions)
    combinations = t(combn(opinions, 2))
    distances = apply(combinations, 1, distance)
    return(table(distances))
}

distances_from_main_peak = function(final_opinions){
    main_peak = names(final_opinions[which(as.numeric(final_opinions) == max(final_opinions))])
    opinions = names(final_opinions)
    other_peaks = opinions[!opinions == main_peak]
    combinations_main_peak = matrix(c(other_peaks, rep(main_peak, length(other_peaks))), ncol=2)
    main_peak_distances = apply(combinations_main_peak, 1, distance)
    return(table(main_peak_distances))
}

all_distances = function(matrix_final_opinions){
    all_final_opinions = apply(unname(matrix_final_opinions), 1, paste, collapse="")
    combinations = t(combn(all_final_opinions, 2))
    distances = apply(combinations, 1, distance)
    return(table(factor(distances, levels=0:ncol(matrix_final_opinions))))
}
