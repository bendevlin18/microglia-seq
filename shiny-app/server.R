
######################################
###                                ###
###     server.R file for          ###
###       microglia-seq            ###
###        application             ###
###                                ###
######################################



server <- function (input, output, session) {
  
  updateSelectizeInput(session, 'entered_genes', choices = genes, server = TRUE)
  
  output$welcome_info <- 
    renderUI({
      c('welcome to the website!')
      
      
    })
  
  output$landing_text <- renderUI({
    
    pub_url <- a('Generation of a microglial developmental index 
                 in mice and in humans reveals a sex difference in 
                 maturation and immune reactivity', href = 'https://onlinelibrary.wiley.com/doi/abs/10.1002/glia.23176', target = '_blank')
    
    x <- list(tags$h1('Welcome to Microglia-Seq!'),
              tags$h3('This website contains interactive visualizations of the data presented in
                       the original publication:'),
              tags$h4(pub_url),
              tags$p('The tabs to the left contain interactive visualizations of the microglial developmental index (MDI) presented in the paper. As
                      well as various plots displaying RPKM values of a list of searched genes (top left search box).'),
              tags$p('All output tables can be downloaded as .csv files for additional analysis!'))
    tagList(x)

    
    
  })
  
  
  output$mdi_plot <- renderPlotly({
    
    df3$age <- factor(df3$age, levels = c('E18', 'P4', 'P14', 'P60'))
    df3$sex <- factor(df3$sex, levels = c('F', 'M'))
    df3$tx <- factor(df3$tx, levels = c('SAL', 'LPS'))
    
    df_grouped <- df3 %>%
      group_by(age, sex, tx) %>%
      summarise(sem = std_err(index), index = mean(index))
    
    p <- ggplot(df_grouped, aes(x = interaction(age,tx), y = index, color = sex, shape = sex), fill = interaction(sex, age))+
      geom_point(data = df3, aes(x = interaction(age,tx), y = index, color = sex), size = 4, position = position_dodge(.5))+
      geom_errorbar(aes(ymin = index - sem, ymax = index + sem), color = 'black', width = .2, size = 1.2, position=position_dodge(.5))+
      graph_theme_settings+
      xlab('\nAge/Treatment')+
      ylab('MDI')+
      ggtitle('Microglial Developmental Index')
    
   if (input$p_vals_mdi){p <- p + stat_compare_means(data = df3, method = 't.test', aes(group = sex), label = 'p.format')}
    
    fig <- ggplotly(p)
    fig
    
  })
  
  
  output$mdi_heatmap <- renderPlotly({
    
    howmany <- input$n_heatmap_genes
    
    ## cleaning up order and factors so that the genes show up in intuitive order in the heatmap
    df4$direction <- factor(df4$direction, levels = c('UP', 'DOWN'))
    df_trunc <- rbind(head(df4, howmany), map_df(tail(df4, howmany), rev))
    df_trunc$abs_val <- abs(df_trunc$valence)
    df_trunc <- map_df(df_trunc[order(df_trunc$abs_val),], rev)
    df_trunc$gene <- factor(df_trunc$gene, levels = rev(df_trunc$gene))
  
    p <- ggplot(df_trunc, aes(y = gene, fill = valence, x = direction))+
      graph_theme_settings+
      xlab('Direction')+
      ylab('Gene')+
      geom_tile(color = 'black')+
      scale_fill_gradient2(midpoint = 0, low = "blue", mid = 'white', high = "red")
    
    fig <- ggplotly(p)
    fig
    
  })
  
  output$mdi_genes <- renderDT(

  df4, options = list(lengthMenu = c(5, 15, 25), pageLength = 15)
  )
  
  output$mdi_gene_download <- downloadHandler(
    
    filename = 'mdi_gene_list_total.csv',
    content = function(file) {
    
      df4
      
      write.csv(df4, file)
      
    })
    
  
  output$summary_table <- renderDT({
    
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(rpkm), rpkm = mean(rpkm), n = n())
    
    
    searched_gene_grouped
    
    })
  
  output$summary_table_download <- downloadHandler(
    
    filename = function() {paste('summary_output', '_', list(input$entered_genes), '.csv', sep = '')},
    content = function(file) {
      goi = input$entered_genes
      
      searched_gene <- df[df$gene %in% goi,]
      
      searched_gene_grouped <- searched_gene %>%
        group_by(gene, age, sex) %>%
        summarise(sem = std_err(rpkm), rpkm = mean(rpkm), n = n())

      write.csv(searched_gene_grouped, file)
      
    })

  
  
  output$bar <- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(rpkm), rpkm = mean(rpkm), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    
    p <- ggplot(searched_gene, aes(age, rpkm, fill = interaction(sex, gene)))+
      geom_bar(stat = 'identity', data = searched_gene_grouped, position = position_dodge(.9))+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = rpkm-sem, ymax = rpkm + sem), width = 0.2, position=position_dodge(.9))+
      facet_wrap(~gene)+
      graph_theme_settings
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = age:sex), label = 'p.format')}
    if (input$ind_points) {p <- p + geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
                                                color = 'black', alpha = 0.7, shape = 21)}
    #if (input$checked_groups){
    # my_comparisons <- list(input$checked_groups)
    # p <- p + stat_compare_means(method = 't.test', comparisons = my_comparison)}
    
    fig <- ggplotly(p)
    fig
  
  })
  
  
  output$violin <- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(rpkm), rpkm = mean(rpkm), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    p <- ggplot(searched_gene, aes(age, rpkm, fill = interaction(sex, gene)))+
      geom_violin(trim = FALSE)+
      facet_wrap(~gene)+
      graph_theme_settings
    
    if (input$p_vals_violin){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    if (input$ind_points_violin) {p <- p + geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
                                                color = 'black', alpha = 0.9)}
    
    fig <- ggplotly(p)
    fig
    
  })
  
  
  output$line<- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(rpkm), rpkm = mean(rpkm), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    p <- ggplot(searched_gene, aes(age, rpkm, group = interaction(sex, gene), linetype = sex))+
      geom_line(size = 1, data = searched_gene_grouped, aes(age, rpkm, group = interaction(sex, gene), color = gene, linetype = sex))+
      geom_point(data = searched_gene_grouped, aes(age, rpkm, group = interaction(sex, gene)), size = 2, color = 'black')+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = rpkm-sem, ymax = rpkm + sem), width = 0.1, color = 'black')+
      graph_theme_settings
    
    if (input$p_vals_line){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    
    fig <- ggplotly(p)
    fig
    
  })
  
  output$dotplot<- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(rpkm), rpkm = mean(rpkm), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    mid <- mean(searched_gene_grouped$rpkm)
    
    p <- ggplot(searched_gene_grouped, aes(x = gene, y = age:sex))+
      geom_point(aes(size = sem, color = rpkm))+
      scale_size(range = c(1,15))+
      scale_color_gradient2(midpoint = mid, low = "blue", mid = 'white', high = "red")+
      coord_flip()+
      dotplot_theme
    
    fig <- ggplotly(p)
    fig
    
  })
  

  
  output$anova_table <- renderDT({
    
    if(input$checked_groups == 'All Ages') {
      goi = input$entered_genes
      aov_output <- aov(df[df$gene == goi[1], ]$rpkm ~ df[df$gene == goi[1], ]$sex * 
                          df[df$gene == goi[1], ]$age)
      aov_df <- data.frame(unclass(summary(aov_output)))
      aov_df['gene'] <- c(goi[1], goi[1], goi[1], goi[1])
      aov_df['labels'] <- c('sex', 'age', 'interaction (sex * age)', 'Residuals')
      row.names(aov_df) <- c('All Ages ', 'All Ages', 'All Ages  ', 'All Ages   ')
      aov_df <- aov_df[, c(7, 1, 2, 3, 4, 5)]
      aov_df <- format(aov_df, digits = 4)
      
    } else {
      
      age_input <- input$checked_groups
      goi = input$entered_genes
      
      df_subset <- df[df$age == age_input[1] | df$age == age_input[2], ]
      
      aov_output <- aov(df_subset[df_subset$gene == goi[1], ]$rpkm ~ df_subset[df_subset$gene == goi[1], ]$sex * 
                          df_subset[df_subset$gene == goi[1], ]$age)
      aov_df <- data.frame(unclass(summary(aov_output)))
      aov_df['gene'] <- c(goi[1], goi[1], goi[1], goi[1])
      aov_df['labels'] <- c('sex', 'age', 'interaction (sex * age)', 'Residuals')
      row.names(aov_df) <- c(paste(age_input[1], '&', age_input[2], ' '), paste(age_input[1], '&', age_input[2]), paste(age_input[1], '&', age_input[2], '  '), paste(age_input[1], '&', age_input[2], '   '))
      aov_df <- aov_df[, c(7, 1, 2, 3, 4, 5)]
      aov_df <- format(aov_df, digits = 4)
    }
    
    aov_df
    
  })
  
  
  output$anova_download <- downloadHandler(
    
    filename = function() {paste('anova_output', '_', input$entered_genes[1], '_', input$checked_groups[1], '_', input$checked_groups[2], '.csv', sep = '')},
    content = function(file) {
    if(input$checked_groups == 'All Ages') {
      goi = input$entered_genes
      aov_output <- aov(df[df$gene == goi[1], ]$rpkm ~ df[df$gene == goi[1], ]$sex * 
                          df[df$gene == goi[1], ]$age)
      aov_df <- data.frame(unclass(summary(aov_output)))
      aov_df['gene'] <- c(goi[1], goi[1], goi[1], goi[1])
      aov_df['labels'] <- c('sex', 'age', 'interaction (sex * age)', 'Residuals')
      row.names(aov_df) <- c('All Ages ', 'All Ages', 'All Ages  ', 'All Ages   ')
      aov_df <- aov_df[, c(7, 1, 2, 3, 4, 5)]
      aov_df <- format(aov_df, digits = 4)
    } else {
      age_input <- input$checked_groups
      goi = input$entered_genes
      
      df_subset <- df[df$age == age_input[1] | df$age == age_input[2], ]
      
      aov_output <- aov(df_subset[df_subset$gene == goi[1], ]$rpkm ~ df_subset[df_subset$gene == goi[1], ]$sex * 
                          df_subset[df_subset$gene == goi[1], ]$age)
      aov_df <- data.frame(unclass(summary(aov_output)))
      aov_df['gene'] <- c(goi[1], goi[1], goi[1], goi[1])
      aov_df['labels'] <- c('sex', 'age', 'interaction (sex * age)', 'Residuals')
      row.names(aov_df) <- c(paste(age_input[1], '&', age_input[2], ' '), paste(age_input[1], '&', age_input[2]), paste(age_input[1], '&', age_input[2], '  '), paste(age_input[1], '&', age_input[2], '   '))
      aov_df <- aov_df[, c(7, 1, 2, 3, 4, 5)]
      aov_df <- format(aov_df, digits = 4)
    }
    aov_df
    write.csv(aov_df, file)
    
  })
  
  
  output$open_target_link <- renderUI({
    
    goi <- input$entered_genes
    url <- a(goi[1], href=df2$ens.id[df2$common.name == goi[1]][1], target = '_blank')
    tagList('Open Targets Link Gene #1: ', url)
    
  })
  
  output$open_target_link2 <- renderUI({
    
    goi <- input$entered_genes
    url2 <- a(goi[2], href=df2$ens.id[df2$common.name == goi[2]][1], target = '_blank')
    tagList('Open Targets Link Gene #2: ', url2)
    
  })
  
  output$open_target_link3 <- renderUI({
    
    goi <- input$entered_genes
    url3 <- a(goi[3], href=df2$ens.id[df2$common.name == goi[3]][1], target = '_blank')
    tagList('Open Targets Link Gene #3: ', url3)
    
  })
  
  output$open_target_link4 <- renderUI({
    
    goi <- input$entered_genes
    url4 <- a(goi[4], href=df2$ens.id[df2$common.name == goi[4]][1], target = '_blank')
    tagList('Open Targets Link Gene #4: ', url4)
    
  })

  
  output$website_info_text <- renderUI({
    
    github_url <- a('Github', href = 'https://github.com/bendevlin18/microglia-seq', target = '_blank')
    
    z <- list(tags$h3('All of the code used to analyze the data and produce this website can be found here: ', github_url),
              tags$p('If you run into any issues or have any questions about the website or the data contained within, contact us!'),
              tags$p('staci.bilbo@duke.edu'),
              tags$p('benjamin.devlin@duke.edu'))
    tagList(z)
    
    
    
  })
  
  
  
  
  
  
  
  
  
}