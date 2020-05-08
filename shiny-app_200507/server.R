
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
    
  
  output$summary_table <- renderDT({
    
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(expression), expression = mean(expression), n = n())
    
    searched_gene_grouped})
  

  
  output$bar <- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(expression), expression = mean(expression), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    
    p <- ggplot(searched_gene, aes(age, expression, fill = interaction(sex, gene)))+
      geom_bar(stat = 'identity', data = searched_gene_grouped, position = position_dodge(.9))+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = expression-sem, ymax = expression + sem), width = 0.2, position=position_dodge(.9))+
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
      summarise(sem = std_err(expression), expression = mean(expression), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    p <- ggplot(searched_gene, aes(age, expression, fill = interaction(sex, gene)))+
      geom_violin(trim = FALSE)+
      facet_wrap(~gene)+
      graph_theme_settings
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    if (input$ind_points) {p <- p + geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
                                                color = 'black', alpha = 0.9)}
    
    fig <- ggplotly(p)
    fig
    
  })
  
  
  output$line<- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(expression), expression = mean(expression), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    p <- ggplot(searched_gene, aes(age, expression, group = interaction(sex, gene), linetype = sex))+
      geom_line(size = 1, data = searched_gene_grouped, aes(age, expression, group = interaction(sex, gene), color = gene, linetype = sex))+
      geom_point(data = searched_gene_grouped, aes(age, expression, group = interaction(sex, gene)), size = 2, color = 'black')+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = expression-sem, ymax = expression + sem), width = 0.1, color = 'black')+
      graph_theme_settings
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    
    fig <- ggplotly(p)
    fig
    
  })
  
  output$dotplot<- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df[df$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(expression), expression = mean(expression), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    mid <- mean(searched_gene_grouped$expression)
    
    p <- ggplot(searched_gene_grouped, aes(x = gene, y = age:sex))+
      geom_point(aes(size = sem, color = expression))+
      scale_size(range = c(1,15))+
      scale_color_gradient2(midpoint = mid, low = "blue", mid = 'white', high = "red")+
      coord_flip()+
      dotplot_theme
    
    fig <- ggplotly(p)
    fig
    
  })
  

  
  output$anova_table <- renderDT({
    
    if(input$twoway_anova) {
      goi = input$entered_genes
      aov_output <- aov(df[df$gene == goi[1], ]$expression ~ df[df$gene == goi[1], ]$sex + 
                          df[df$gene == goi[1], ]$age)
      aov_df <- data.frame(unclass(summary(aov_output)))
      aov_df['labels'] <- c('sex', 'age', 'Residuals')
      aov_df['gene'] <- c(goi[1], goi[1], goi[1])
      aov_df <- aov_df[, c(7, 6, 1, 2, 3, 4, 5)]
      
    } else {
      
      age_input <- input$checked_groups
      goi = input$entered_genes
      
      df_subset <- df[df$age == age_input[1] | df$age == age_input[2], ]
      
      aov_output <- aov(df_subset[df_subset$gene == goi[1], ]$expression ~ df_subset[df_subset$gene == goi[1], ]$sex + 
                          df_subset[df_subset$gene == goi[1], ]$age)
      aov_df <- data.frame(unclass(summary(aov_output)))
      aov_df['labels'] <- c('sex', 'age', 'Residuals')
      aov_df['gene'] <- c(goi[1], goi[1], goi[1])
      aov_df['age'] <- c(age_input[1], age_input[2], '')
      aov_df <- aov_df[, c(8, 7, 6, 1, 2, 3, 4, 5)]
    }
    
    aov_df
    
  })
  
  
  output$open_target_link <- renderUI({
    
    goi <- input$entered_genes
    
    url <- a(goi[1], href=df2$ens.id[df2$common.name == goi[1]][1], target = '_blank')
    url2 <- a(goi[2], href=df2$ens.id[df2$common.name == goi[2]][1], target = '_blank')
    url3 <- a(goi[3], href=df2$ens.id[df2$common.name == goi[3]][1], target = '_blank')
    url4 <- a(goi[4], href=df2$ens.id[df2$common.name == goi[4]][1], target = '_blank')
    tagList('Open Targets Links: ', list(url, url2, url3, url4))
    
  })
  
}