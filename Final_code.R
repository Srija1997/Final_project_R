library(tidyverse)
library(shiny)
library(ggplot2)
library(colourpicker)
library(DT)
library(ggbeeswarm)
library(RColorBrewer)
library(pheatmap)

options(shiny.maxRequestSize = 30*1024^2)
options(warn=-1)



ui <- fluidPage(
  headerPanel("Srija Chillamcherla: Final Project"),
  tabsetPanel(
    tabPanel('Samples',
             sidebarLayout(
               sidebarPanel(
                 fileInput("exp_file", label="Load samples metadata here:", accept = c(".csv")),
                 p ('A summary table and plots for the samples from the study.'),
                 submitButton("Submit", icon("check")),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            p ('A summary of the metdata of the samples'),
                            tableOutput('summary')),
                   tabPanel("Table",
                            p ('A table with sortable columns'),
                            DT::dataTableOutput('table')),
                   tabPanel("Plots",
                            p ('A histogram plot can be generated with from the chosen column.'),
                            sidebarPanel(
                              selectInput("plot_col", "Select the column for the plot", 'Age_of_death', choices = c('Age_of_death', 'PMI', 'RIN', 'Seq_reads')) 
                            ),
                            submitButton("Submit"),
                            plotOutput('plot'))
                 )
               )
             )  
    ),
    tabPanel('Normalised counts matrix',
             sidebarLayout(
               sidebarPanel(
                 fileInput("norm_file", label="Load normalised counts matrix:", accept = c(".csv")),
                 p ('Exploring and visualizing counts matrices can aid in selecting gene count filtering strategies and understanding counts data structure.'),
                 sliderInput(inputId = "variance", min = 0, max = 1,
                             label = "Select the variance level to select the genes:", value=0.5, step = 0.01),
                 sliderInput(inputId = "samples_nonzero", min = 1, max = 69,
                             label = "Select the number of genes with non-zero samples:", value=35, step = 2),
                 submitButton("Submit", icon("check")),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            p ('Table summarising the genes after the filtering'),
                            tableOutput('filtered_summary')),
                   tabPanel("Scatter plots",
                            p ('Diagnostic scatter plots for the median count vs variance, median count vs number of zeros of the genes filtered (darker colour), omitted genes (lighter colour)'),
                            fluidRow(
                              splitLayout(cellWidths = c("50%", "50%"),plotOutput("varplot"),plotOutput("naplot"))
                            )),
                   tabPanel("Heatmap",
                            p ('Clustered heatmap of counts remaining after filtering'),
                            plotOutput('filtered_heatmap')),
                   tabPanel("Principal Component Analysis",
                            p ('PCA plot'),
                            sidebarPanel(
                              selectInput("x_pc", "Select the principal component for X-axis", '1', choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                              selectInput("y_pc", "Select the principal component for Y-axis", '2', choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                              submitButton("Submit"),
                            ),
                            plotOutput('filtered_pca'))
                 )
               )
             )
    ),
    tabPanel('Differential expression results',
             sidebarLayout(
               sidebarPanel(
                 fileInput("diff_file", label="Load differential expression results:", accept = c(".txt")),
                 p ('Differential expression identifies which genes, if any, are implicated in a specific biological comparison.'),
                 submitButton("Submit", icon("check")),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Results",
                            p ('Sortable table displaying differential expression results'),
                            DT::dataTableOutput('diff_results')),
                   tabPanel("Plot",
                            sidebarPanel(
                              p ('A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.'),
                              radioButtons("x_axis", "Choose the column for X axis:", "log2FoldChange",
                                           choices = c('baseMean',
                                                       'log2FoldChange',
                                                       'stat',
                                                       'pvalue',
                                                       'padj',
                                                       'lfcSE')),
                              radioButtons("y_axis", "Choose the column for Y axis:", "pvalue",
                                           choices = c('baseMean',
                                                       'log2FoldChange',
                                                       'stat',
                                                       'pvalue',
                                                       'padj',
                                                       'lfcSE')),
                              colourInput("base_point_col", 'Base point color', "#1A6191"),
                              colourInput("highlight_col", 'Highlight point color', "#9223A1"),
                              sliderInput(inputId = "p_magnitude", min = -12, max = 10,
                                          label = "Select the magnitude of the p adjusted coloring:", value = 0, step = 0.5),
                              submitButton("Submit", icon("check"))
                            ),
                            
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table",
                                         tableOutput('diff_table')),
                                tabPanel("Plot",
                                         plotOutput('diff_plot'))
                              )
                            ),
                   ),
                 )
               )
             )
    ),
    tabPanel('Individual gene expression',
             tabPanel('Samples',
                      sidebarLayout(
                        sidebarPanel(
                          p ('Summary and visualizing individual gene counts is sometimes useful for examining or verifying patterns identified by differential expression analysis.'),
                          fileInput("gene_file", label="Load samples file here:", accept = c(".csv")),
                          fileInput("norm_file2", label="Load normalised counts matrix:", accept = c(".csv")),
                          textInput(inputId = "search_term",label = "Search gene here:", placeholder = "ENSG00000000419.8"),
                          submitButton("Submit", icon("check")),
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Visualisation",
                                     p ('The gene can be visualised in the form of a bar plot, boxplot, violin plot or beeswarm plot.'),
                                     sidebarPanel(
                                       selectInput(inputId = "gene_col", label = "Select the column for the plot", 'Age_of_death', choices = c('Age_of_death', 'PMI', 'RIN', 'Seq_reads')),
                                       selectInput(inputId = "plot_type", label = "Select the type of the plot", 'Boxplot', choices = c('Bar plot', 'Boxplot', 'Violin plot', 'Beeswarm plot')),
                                       submitButton("Submit")
                                     ),
                                     plotOutput('gene_plot'))
                          )
                        )
                      )  
             ),),
  )
)


server <- function(input, output, session) {
  load_data <- reactive({
    file <- input$exp_file
    if(is.null(file)){
      return(NULL)
    }
    return(read.csv(file$datapath))
  })
  
  draw_table <- function(df){
    file <- input$exp_file
    if(is.null(file)){return(NULL)}
    RM<-as.matrix(colMeans(df[,9:12])) 
    colnames(RM) <- "Mean"
    return(RM)
  }
  
  plot_hist <- function(meta,col){
    hist <- ggplot(meta,aes(x = get(col)))+geom_histogram(bins=10, color = "darkblue", fill = "lightblue", alpha = 0.8)+labs(title='Histogram of the selected column')+theme_minimal()+theme_classic()+theme(legend.position = "bottom")+xlab(col)+ylab('Count')
    hist <- hist + theme(text=element_text(size=16))
    return(hist)
  }
  
  norm_data <- reactive({
    options(shiny.maxRequestSize = 30*1024^2)
    norm <- input$norm_file
    if(is.null(norm)){return(NULL)}
    return(read.csv(norm$data))
  })
  
  norm_summary <- function(dat, var, zero){
    norm <- input$norm_file
    if(is.null(norm)){return(NULL)}
    dat$row_var = apply(dat[-1], 1, var)
    dat$count_na <- apply(dat[-1] == 0, 1, sum)
    f <- dplyr::filter(dat, dat$row_var >= quantile(row_var, var) & dat$count_na >= zero)
    unfiltered_genes <- nrow(dat)-nrow(f)
    df <- data.frame(Description = c("Total number of genes", "Total number of samples", "Number of filtered genes", "Percentage of filtered genes", "Number of unfiltered genes", "Percentage of unfiltered genes"),
                     Number = c(nrow(dat), ncol(dat), nrow(f), nrow(f)/nrow(dat)*100, unfiltered_genes, unfiltered_genes/nrow(dat)*100)) %>% as_tibble()
    return(df)
  }
  
  
  norm_var_plot <- function(data, var, zero){
    norm <- input$norm_file
    if(is.null(norm)){return(NULL)}
    data$row_median <- apply(data[-1],1, median)
    data$row_var = apply(data[-1], 1, var)
    data$count_na <- apply(data[-1] == 0, 1, sum)
    data$filter = "unfiltered"
    data$filter[data$row_var >= var & data$count_na >= zero] = "filtered"
    plot_var <- ggplot()+
      geom_point(data, mapping = aes(log10(row_median), (row_var), color=filter))
    plot_var <- plot_var + scale_y_log10()
    return(plot_var)
    
  }
  
  norm_na_plot <- function(data, var, zero){
    norm <- input$norm_file
    if(is.null(norm)){return(NULL)}
    data$row_median <- apply(data[-1],1, median)
    data$row_var = apply(data[-1], 1, var)
    data$count_na <- apply(data[-1] == 0, 1, sum)
    data$filter = "unfiltered"
    data$filter[data$row_var >= var & data$count_na >= zero] = "filtered"
    plot_na <- ggplot()+
      geom_point(data, mapping = aes(log10(row_median), (count_na), color=filter))
    return(plot_na)
    
  }
  
  norm_heatmap <- function(data, var, zero){
    norm <- input$norm_file
    if(is.null(norm)){return(NULL)}
    data$row_median <- apply(data[-1],1, median)
    data$row_var = apply(data[-1], 1, var)
    data$count_na <- apply(data[-1] == 0, 1, sum)
    data$filter = "unfiltered"
    data$filter[data$row_var >= var & data$count_na >= zero] = "filtered"
    filtered_genes <- subset(data[ ,!(colnames(data) %in% c("row_var","row_median","count_na"))], filter=="filtered")
    filtered_genes <- filtered_genes[ ,!(colnames(filtered_genes) %in% c("filter"))]
    
    p <- pheatmap(as.matrix(filtered_genes[-1]), scale = 'row', show_rownames = F, col = brewer.pal(11, "Set3"))
    return(p)
    
  }
  
  filtered_pca_plot <- function(pca_data, pc1, pc2){
    norm <- input$norm_file
    if(is.null(norm)){return(NULL)}
    pca_results <- prcomp(scale(t(pca_data[-1])), center=FALSE, scale=FALSE)
    percent_var <- pca_results$sdev^2 / sum( pca_results$sdev^2 )
    plot <- as_tibble(pca_results$x) %>%
      ggplot(aes(x=get(pc1), y = get(pc2)))+
      geom_point()+
      xlab(sprintf("%s (Variance: %s)", input$x_pc, percent_var[strtoi(substr(input$x_pc, 3, 3))]))+
      ylab(sprintf("%s (Variance: %s)", input$y_pc, percent_var[strtoi(substr(input$y_pc, 3, 3))]))+
      ggtitle("Principal Component Analysis of the chosen PCs")
    return(plot)
  }
  
  diff_data <- reactive({
    options(shiny.maxRequestSize = 80*1024^2)
    diff <- input$diff_file
    if(is.null(diff)){return(NULL)}
    return(read.csv(diff$datapath, header = TRUE))
  })
  
  diff_results_table <- function(data){
    diff <- input$diff_file
    if(is.null(diff)){return(NULL)}
    data <- dplyr::filter(data, data$pvalue < 0.05)
    data$pvalue <- lapply(data$pvalue, formatC, digits = 9)
    return(data)
  }
  
  difftable <- function(df, slider) {
    diff <- input$diff_file
    if(is.null(diff)){return(NULL)}
    df <- dplyr::filter(df, df$padj < (1*(10**slider)))
    df$pvalue <- lapply(df$pvalue, formatC, digits = 9)
    df$padj <- lapply(df$padj, formatC, digits = 9)
    
    return(df)
  }
  
  volcano_plot <- function(df, x_col, y_col, slider, color1, color2) {
    diff <- input$diff_file
    if(is.null(diff)){return(NULL)}
    y_lab = sprintf("-log10(%s)", y_col)
    col_label = sprintf("padj < (1*(10^%s))", slider)
    plot <-  ggplot(df, aes(x = get(x_col), y = -log10(get(y_col)), color = padj < (1*(10**(slider))))) +
      geom_point() +
      scale_colour_manual(name = col_label, values = setNames(c(color1, color2, 'grey'),c(T, F, NA))) +
      xlab(x_col) +
      ylab(y_col) +
      theme_linedraw()
    return(plot)
  }
  
  sample_data <- reactive({
    options(shiny.maxRequestSize = 80*1024^2)
    genefile <- input$gene_file
    if(is.null(genefile)){return(NULL)}
    return(read.csv(genefile$datapath, header = TRUE))
  })
  
  norm_data2 <- reactive({
    options(shiny.maxRequestSize = 30*1024^2)
    norm2 <- input$norm_file2
    if(is.null(norm2)){return(NULL)}
    return(read.csv(norm2$data))
  })
  
  geneplot <- function(sample_dat, norm_dat2, col, term, gene_plot_type){
    genefile <- input$gene_file
    if(is.null(genefile)){return(NULL)}
    
    norm2 <- input$norm_file2
    if(is.null(norm2)){return(NULL)}
    
    norm_dat2 <- norm_dat2 %>% remove_rownames %>% column_to_rownames(var="X")
    
    sample_dat$Age_of_death   <- cut(sample_dat$Age_of_death,seq(40,100,10),right=FALSE,labels=c("Age group:40-49","Age group:50-59","Age group:60-69","Age group:70-79","Age group:80-89","Age group:90-99"))
    sample_dat$PMI <- cut(sample_dat$PMI,seq(0,35,7),right=FALSE,labels=c("PMI:0-6","PMI:7-14","PMI:15-21","PMI:22-28","PMI:29-35"))
    sample_dat$RIN <- cut(sample_dat$RIN, seq(6,10,1), right=FALSE,labels=c("RIN:6-6.9","RIN:7-7.9","RIN:8-8.9","RIN:9-10"))
    sample_dat$Seq_reads <- cut(sample_dat$Seq_reads, seq(38420004,167044880,20000000), right=FALSE, labels=c("Group1", "Group2", "Group3", "Group4", "Group5","Group6"))
    
    col_variable <- input$gene_col
    gene <- input$search_term
    
    xcol <- unlist(select(sample_dat, col_variable))
    ycol <- unlist(norm_dat2[gene,]) 
    
    g <- ggplot(mapping=aes(x=xcol,y=ycol))+xlab(col_variable)+ylab(gene)
    
    if(gene_plot_type == "Bar plot"){
      barplot <- g + geom_bar(stat = 'identity', fill = '#69b3a2')
      return(barplot)
    }
    
    if(gene_plot_type == "Boxplot"){
      boxplot <- g + geom_boxplot(fill = 'pink', color = 'steelblue')
      return(boxplot)
    }
    
    
    if(gene_plot_type == "Violin plot"){
      violinplot <- g + geom_violin(fill="#E69F00", col = "#999999")
      return(violinplot)
    }
    
    
    if(gene_plot_type == "Beeswarm plot"){
      beeswarmplot <- g + geom_beeswarm(col = "red", pch = 18,cex = 1.5, method = "hex")
      return(beeswarmplot)
    }
    
  }
  
  #Tab 1
  output$summary <- renderTable(draw_table(load_data()), rownames = TRUE)
  output$table <- DT::renderDataTable((load_data()))
  output$plot <- renderPlot(plot_hist(load_data(),input$plot_col))
  
  # Tab 2
  output$filtered_summary <- renderTable(norm_summary(norm_data(), input$variance, input$samples_nonzero))
  output$varplot <- renderPlot(norm_var_plot(norm_data(), input$variance, input$samples_nonzero))
  output$naplot <- renderPlot(norm_na_plot(norm_data(), input$variance, input$samples_nonzero))
  output$filtered_heatmap <- renderPlot(norm_heatmap(norm_data(), input$variance, input$samples_nonzero),  height = 600, width = 600)
  output$filtered_pca <- renderPlot(filtered_pca_plot(norm_data(), input$x_pc, input$y_pc))
  
  #Tab 3
  output$diff_results <- DT::renderDataTable((diff_data()))
  output$diff_table <- renderTable(difftable(diff_data(), input$p_magnitude))
  output$diff_plot <- renderPlot(volcano_plot(diff_data(), input$x_axis, input$y_axis, input$p_magnitude, input$base_point_col, input$highlight_col), height = 600, width = 600)
  
  #Tab 4
  output$gene_plot <- renderPlot(geneplot(sample_data(), norm_data2(), input$gene_col, input$search_term, input$plot_type))
}

shinyApp(ui, server)
