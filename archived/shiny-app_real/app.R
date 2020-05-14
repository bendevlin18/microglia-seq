
##### Load in necessary packages #####

library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(plotly)
library(shinythemes)

### Load in the data and generate objects ###

df <- as_tibble(read.csv('GSE99622_hanamsagar2017_cleaned_melted.csv'))
df2 <- as_tibble(read.csv('gene_ensembl_ids_opentarget.csv'))

genes <- unique(df$gene)
std_err <- function(x) sd(x) / sqrt(length(x))

x <- unique(df$age)


ui <- fluidPage(theme = shinytheme('journal'),
                
                titlePanel('Microglia Development RNASeq'),
                
                sidebarLayout(
                  sidebarPanel(
                    selectizeInput('entered_genes', label = 'Search for one or multiple genes', 
                                   choices = genes, multiple = TRUE, options = list(maxItems = 5, 
                                   placeholder = 'Example: Cx3cr1')),
                      checkboxGroupInput('checked_groups', label = 'Choose ages to compare', choices = x),
                      checkboxInput('p_vals', label = 'Turn P values ON (For Bar, Violin, and single gene Line plots ONLY)'),
                      checkboxInput('ind_points', label = 'Turn Individual Points ON'),
                      htmlOutput('open_target_link')),
                  
                  mainPanel(
                    tabsetPanel(
                      tabPanel('Bar Plot', plotlyOutput('bar', width = 'auto', height = 900)),
                      tabPanel('Violin Plot', plotlyOutput('violin', width = 'auto', height = 900)),
                      tabPanel('Dot Plot', plotlyOutput('dotplot', width = 'auto', height = 900)),
                      tabPanel('Line Plot', plotlyOutput('line', width = 'auto', height = 900)),
                      tabPanel('Summary Data Table', tableOutput('summary_table'))
                      
                    ))))



server <- function (input, output) {
  
  
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
    
    p <- ggplot(searched_gene, aes(age, expression, fill = sex:gene))+
             geom_bar(stat = 'identity', data = searched_gene_grouped, position = position_dodge(.9))+
             geom_errorbar(data = searched_gene_grouped, aes(ymin = expression-sem, ymax = expression + sem), width = 0.2, position=position_dodge(.9))+
             xlab('\r\n Age/Treatment')+
             ylab('Average Expression')+
             facet_wrap(~gene)+
             theme_classic()+
             theme(rect = element_rect(fill = 'transparent'), 
                   text = element_text(size = 20), 
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
    
    p <- ggplot(searched_gene, aes(age, expression, fill = sex:gene))+
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
    geom_line(size = 1, data = searched_gene_grouped, aes(age, expression, group = sex:gene, color = gene, linetype = sex))+
    geom_point(data = searched_gene_grouped, aes(age, expression, group = sex:gene), size = 2, color = 'black')+
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

output$open_target_link <- renderUI({
  
  goi <- input$entered_genes
  
  url <- a(goi[1], href=df2$ens.id[df2$common.name == goi[1]][1], target = '_blank')
  
  tagList('Open Targets Link: ', url)
  
})

}
  
  shinyApp(ui = ui, server = server)







