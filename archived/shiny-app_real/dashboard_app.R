
##### Load in necessary packages #####

library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(plotly)
library(shinythemes)



ui <- dashboardPage(skin = 'black',
                
                dashboardHeader(title = 'Microglia Development RNASeq', titleWidth = 375),
                
                dashboardSidebar(
                    selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for one or multiple genes',
                                   options = list(maxItems = 5, placeholder = 'Example: Cx3cr1')),
                    sidebarMenu(
                      menuItem('Plots', tabName = 'Plots', icon = icon('bar-chart-o'), startExpanded = FALSE,
                               menuSubItem('Bar Plot', tabName = 'BarPlot'),
                               menuSubItem('Violin Plot', tabName = 'ViolinPlot'),
                               menuSubItem('Line Plot', tabName = 'LinePlot'),
                               menuSubItem('Dot Plot', tabName = 'DotPlot')),
                      menuItem('Info', tabName = "SummaryDataTable", icon = icon('info')),
                      menuItem('Stats', tabName = 'ANOVATable', icon = icon('th'))
                    ),
                    
                    checkboxGroupInput('checked_groups', label = 'Choose two ages to compare (Two-Way Anova)', 
                                       choices = c('E18', 'P4', 'P14', 'P60', 'P60 + LPS')),
                    checkboxInput('twoway_anova', label = 'Compare Across All Ages (Two-Way Anova)'),
                    checkboxInput('p_vals', label = 'Turn P values ON (For Bar, Violin, and single gene Line plots ONLY)'),
                    checkboxInput('ind_points', label = 'Turn Individual Points ON'),
                    htmlOutput('open_target_link'),
                    htmlOutput('open_target_link2')
                ),
                
                dashboardBody(
                  tabItems(
                    tabItem(tabName = 'BarPlot', plotlyOutput('bar', width = 'auto', height = 'auto')),
                    tabItem(tabName = 'ViolinPlot', plotlyOutput('violin', width = 'auto', height = 'auto')),
                    tabItem('DotPlot', plotlyOutput('dotplot', width = 'auto', height = 'auto')),
                    tabItem('LinePlot', plotlyOutput('line', width = 'auto', height = 'auto')),
                    tabItem('SummaryDataTable', tableOutput('summary_table')),
                    tabItem(tabName = 'ANOVATable', tableOutput('anova_table'))))
                )


server <- function (input, output, session) {
  
  ### Load in the data and generate objects ###
  df <- as_tibble(read.csv('GSE99622_hanamsagar2017_cleaned_melted.csv'))
  df2 <- as_tibble(read.csv('gene_ensembl_ids_opentarget.csv'))
  
  genes <- unique(df$gene)
  std_err <- function(x) sd(x) / sqrt(length(x))
  
  updateSelectInput(session, 'entered_genes', choices = genes)
  
  
  
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
             xlab('\r\n Age/Treatment')+
             ylab('Average Expression')+
             facet_wrap(~gene)+
             theme_classic()+
             theme(rect = element_rect(fill = 'transparent'), 
                   text = element_text(size = 10), 
                   axis.line = element_line(size = 1),
                   axis.ticks = element_line(size = .6),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   legend.text = element_text(size = 15),
                   axis.ticks.length = unit(.2, 'cm'),
                   panel.background = element_rect(fill = "transparent"),
                   plot.background = element_rect(fill = "transparent", color = NA),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.background = element_rect(color = NA, col = 0),
                   axis.title.y = element_text(margin = margin(t = 0, r = 75, b = 0, l = 0)),
                   axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
                   strip.background = element_rect(color="transparent", fill="transparent", size=5, linetype="solid"),
                   panel.spacing = unit(1, 'cm'))
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = age:sex), label = 'p.format')}
    if (input$ind_points) {p <- p + geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
                                                             color = 'black', alpha = 0.7, shape = 21)}
    #if (input$checked_groups){
     # my_comparisons <- list(input$checked_groups)
     # p <- p + stat_compare_means(method = 't.test', comparisons = my_comparison)}
    
    fig <- ggplotly(p)
    print(fig)
    
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
              xlab('Age/Treatment')+
              ylab('Average Expression')+
              facet_wrap(~gene)+
              theme_classic()+
              theme(rect = element_rect(fill = 'transparent'), 
                    text = element_text(size = 30), 
                    axis.line = element_line(size = 1),
                    axis.ticks = element_line(size = .6),
                    axis.ticks.length = unit(.2, 'cm'),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    panel.background = element_rect(fill = "transparent"),
                    plot.background = element_rect(fill = "transparent", color = NA),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    legend.background = element_rect(color = NA, col = 0),
                    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
                    axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
                    strip.background = element_rect(color="transparent", fill="transparent", size=5, linetype="solid"),
                    panel.spacing = unit(1, 'cm'))
    
    if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
    if (input$ind_points) {p <- p + geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
                                                color = 'black', alpha = 0.9)}
    
    fig <- ggplotly(p)
    print(fig)
    
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
  
  p <- ggplot(searched_gene, aes(age, expression, group = sex:gene, linetype = sex))+
    geom_line(size = 1, data = searched_gene_grouped, aes(age, expression, group = interaction(sex, gene), color = gene, linetype = sex))+
    geom_point(data = searched_gene_grouped, aes(age, expression, group = interaction(sex, gene)), size = 2, color = 'black')+
    geom_errorbar(data = searched_gene_grouped, aes(ymin = expression-sem, ymax = expression + sem), width = 0.1, color = 'black')+
    xlab('Age/Treatment')+
    ylab('Average Expression')+
    theme_classic()+
    theme(rect = element_rect(fill = 'transparent'), 
          text = element_text(size = 30), 
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = .6),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.length = unit(.2, 'cm'),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(color = NA, col = 0),
          legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          strip.background = element_rect(color="transparent", fill="transparent", size=5, linetype="solid"),
          panel.spacing = unit(1, 'cm'))
  
  if (input$p_vals){p <- p + stat_compare_means(method = 't.test', aes(group = sex), label = 'p.format')}
  
  fig <- ggplotly(p)
  print(fig)
  
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
             theme_bw()+
             scale_size(range = c(1,15))+
             scale_color_gradient2(midpoint = mid, low = "blue", mid = 'white', high = "red")+
             coord_flip()+
             theme(
               text = element_text(face = 'bold', size = 20),
               axis.text.x = element_text(angle = 45, hjust = 1),
               axis.title.x = element_blank())
  
  fig <- ggplotly(p)
  print(fig)
  
})

output$summary_table <- renderTable({
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
  
  searched_gene_grouped})

output$anova_table <- renderTable({
  
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
  
  tagList('Open Targets Link: ', url)
  
})

}
  
  shinyApp(ui = ui, server = server)







