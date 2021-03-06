

######################################
###                                ###
###     server.R file for          ###
###       microglia-seq            ###
###           mobile               ###
###        application             ###
###                                ###
######################################



server <- function (input, output, session) {
  
  updateSelectizeInput(session, 'entered_genes', choices = genes, server = TRUE)
  
  output$landing_text <- renderUI({
    
    pub_url <- a('Generation of a microglial developmental index 
                 in mice and in humans reveals a sex difference in 
                 maturation and immune reactivity', href = 'https://onlinelibrary.wiley.com/doi/abs/10.1002/glia.23176')
    
    x <- list(tags$h1('Welcome to Microglia-Seq Mobile!'),
              tags$h3('This website contains interactive visualizations of the data presented in
                       the original publication:'),
              tags$h4(pub_url),
              fluidRow(
                column(width = 12,
                       div(style = "height:35px;background-color: transparent;"))),
              tags$p('The tabs below contain interactive visualizations of the microglial developmental index (MDI) presented in the paper. As
                      well as various plots/tables of TPM values from list of searched genes'),
              fluidRow(
                column(width = 12,
                       div(style = "height:35px;background-color: transparent;"))),
              tags$b('Video Abstract'))
    tagList(x)
    
  })
  
  
  output$mdi_plot <- renderPlotly({
    
    df3_m$age <- factor(df3_m$age, levels = c('E18', 'P4', 'P14', 'P60'))
    df3_m$sex <- factor(df3_m$sex, levels = c('F', 'M'))
    df3_m$tx <- factor(df3_m$tx, levels = c('SAL', 'LPS'))
    
    df_grouped <- df3_m %>%
      group_by(age, sex, tx) %>%
      summarise(sem = std_err(index), index = mean(index))
    
    p <- ggplot(df_grouped, aes(x = interaction(age,tx), y = index, color = sex, shape = sex), fill = interaction(sex, age))+
      geom_point(data = df3_m, aes(x = interaction(age,tx), y = index, color = sex), size = 2, position = position_dodge(.5))+
      geom_errorbar(aes(ymin = index - sem, ymax = index + sem), color = 'black', width = .2, size = 1.2, position=position_dodge(.5))+
      graph_theme_settings+
      xlab('\nAge/Treatment')+
      ylab('MDI')+
      ggtitle('Microglial Developmental Index')
    
    if (input$p_vals){p <- p + stat_compare_means(data = df3_m, method = 't.test', aes(group = sex), label = 'p.signif')}
    
    ggplotly(p) %>%
      layout(legend = list(x = 0.9, y = 1)) %>%
      config(displayModeBar = FALSE)
    
  })
  
  output$mdi_genes <- renderDT(
    
    datatable(df4_m, options = list(sDom  = '<"top"><"bottom">ip'))
  )
  
  
  output$summary_table <- renderDT({
    
    goi = input$entered_genes
    
    searched_gene <- df_m[df_m$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(TPM), TPM = mean(TPM), n = n())
    
    datatable(searched_gene_grouped, options = list(sDom  = '<"top"><"bottom">ip'))
    
  })
  
  output$line_summary_table <- renderDT({
    
    goi = input$entered_genes
    
    searched_gene <- df_m[df_m$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(TPM), TPM = mean(TPM), n = n())
    
    datatable(searched_gene_grouped, options = list(sDom  = '<"top"><"bottom">ip'))
    
  })
  
  output$dot_summary_table <- renderDT({
    
    goi = input$entered_genes
    
    searched_gene <- df_m[df_m$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(TPM), TPM = mean(TPM), n = n())
    
    datatable(searched_gene_grouped, options = list(sDom  = '<"top"><"bottom">ip'))
    
  })
  
  
  output$bar <- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df_m[df_m$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(TPM), TPM = mean(TPM), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    
    p <- ggplot(searched_gene, aes(age, TPM, fill = interaction(sex, gene)))+
      geom_bar(stat = 'identity', data = searched_gene_grouped, position = position_dodge(.9))+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = TPM-sem, ymax = TPM + sem), width = 0.2, position=position_dodge(.9))+
      ggtitle(goi, 'Expression')+
      graph_theme_settings
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.signif')}
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = age:sex), label = 'p.signif')}
    if (input$ind_points) {p <- p + geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
                                                color = 'black', alpha = 0.7, shape = 21)}
    
    ggplotly(p) %>%
      layout(legend = list(x = 0.9, y = 1)) %>%
      config(displayModeBar = FALSE)
    
  })
  
  
  output$line<- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df_m[df_m$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(TPM), TPM = mean(TPM), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    p <- ggplot(searched_gene, aes(age, TPM, group = interaction(sex, gene), linetype = sex))+
      geom_line(size = 1, data = searched_gene_grouped, aes(age, TPM, group = interaction(sex, gene), color = gene, linetype = sex))+
      geom_point(data = searched_gene_grouped, aes(age, TPM, group = interaction(sex, gene)), size = 2, color = 'black')+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = TPM-sem, ymax = TPM + sem), width = 0.1, color = 'black')+
      ggtitle(goi, 'Expression')+
      graph_theme_settings
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.signif')}
    
    ggplotly(p) %>%
      layout(legend = list(x = 0.9, y = 1)) %>%
      config(displayModeBar = FALSE)
    
  })
  
  output$dotplot<- renderPlotly({
    goi = input$entered_genes
    
    searched_gene <- df_m[df_m$gene %in% goi,]
    
    searched_gene_grouped <- searched_gene %>%
      group_by(gene, age, sex) %>%
      summarise(sem = std_err(TPM), TPM = mean(TPM), n = n())
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    mid <- mean(searched_gene_grouped$TPM)
    
    p <- ggplot(searched_gene_grouped, aes(x = gene, y = age:sex))+
      geom_point(aes(size = sem, color = TPM))+
      scale_size(range = c(1,5))+
      scale_color_gradient2(midpoint = mid, low = "blue", mid = 'white', high = "red")+
      coord_flip()+
      dotplot_theme
    
    ggplotly(p) %>%
      layout(legend = list(x = 0.9, y = 1)) %>%
      config(displayModeBar = FALSE)
    
  })
  
  
  output$website_info_text <- renderUI({
    
    z <- list(
      fluidRow(
        column(width = 12,
               div(style = "height:50px;background-color: transparent;"))),
      tags$h3('All of the code used to analyze the data and produce this website can be found here: ', tags$a(
        href = 'https://github.com/bendevlin18/microglia-seq', target = '_blank',
        tags$img(src="github-logo.png", title="GitHub", width="25", height="25"))),
      tags$p('If you run into any issues or have any questions about the website or the data contained within, contact us!'),
      tags$a('Return to Desktop Version', href = 'http://microglia-seq.vm.duke.edu/microglia-seq/shiny-app/', target = '_blank'),
      fluidRow(
        column(width = 12,
               div(style = "height:50px;background-color: transparent;"))),
      fluidRow(
        column(width = 2),
        column(width = 4, 
               tags$img(src="staci.jpg", title="Staci Image", width="225", height="231"),
               tags$p('Staci Bilbo, Ph.D'),
               tags$p('Lab Head'),
               tags$a('staci.bilbo@duke.edu')),
        column(width = 4, 
               tags$img(src="ben.png", title="Ben Image", width="225", height="225"),
               tags$p('Ben Devlin'),
               tags$p('Graduate Student'),
               tags$a('benjamin.devlin@duke.edu')),
        column(width = 2))
      
    )
    
    
    
    tagList(z)
    
    
    
  })
  
  
  
  
  
  
  
  
  
}