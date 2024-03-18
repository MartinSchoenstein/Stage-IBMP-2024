INTERVAL = 200
valuesInterval = c(1, 1 + INTERVAL)


ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(
    width = 250,
    sidebarMenu(id = "sidebar",
     menuItem("PolyA", tabName = "polya", icon = icon("font")),
     menuItem("Uridylation", tabName = "uridyl", icon = icon("magnet")),
     menuItem("Reads", tabName = "reads", icon = icon("layer-group")),
     conditionalPanel(condition = 'input.sidebar == "uridyl"',
        selectInput("gene", "Selectionnez le gène : ", listGene),
        textInput("geneval", "Ou écrivez le :")),
     conditionalPanel(condition = 'input.sidebar == "polya"',
        selectInput("barcode", "Quel barcode ?", c("bc6", "bc9", "bc10", "bc12")))
     
    )),
  dashboardBody(
    tabItems(
      tabItem(tabName = "uridyl",
        fluidRow(
          box(title = "Percentage of Uridylation :", status = "success", solidHeader = TRUE, width = 12, plotOutput("plott", width = "850px"))
  
    )
  ),
      tabItem(tabName = "polya",
        fluidRow(
          box(title = "Distribution of polyA tail lenght : ", status = "success", solidHeader = TRUE, width = 12, plotOutput("plot", width = "650px")),
          
          box(sliderInput("xrange", "Modifiez l'axe des x :", min = 0, max = max(fbindlist$polya_length), value = max(fbindlist$polya_length))),
          
          box(textInput("xrangval", "xmax :", value = max(fbindlist$polya_length)))
        
  )
),

      tabItem(tabName = "reads",
        fluidRow(
         box(title = "Choose a part of the table :", width = 12, conditionalPanel(condition = 'nrow(fbindlist) > 200', sliderInput("intervalReads", "Reads : ", value = valuesInterval, min = 1, max = 5000, width = '100%'))),
         
         box(title = "Distribution of reads length :", status = "success", width = 12, solidHeader = TRUE, collapsible = TRUE, plotOutput("plotreads")),
           
         box(title = "DataTable : ", status = "primary", width = 12, solidHeader = TRUE, collapsible = TRUE, downloadButton("downloadtable", "Download selected table"), p(""), dataTableOutput("table")))
    ) 
  )
))

server = function(input, output, session) {
  reactGDF = reactive({
    geneDataframe = fbindlist%>%mutate(urydil = if_else(add_tail_pct_T>60, "U", "/"))%>%filter(mRNA == input$gene)
  })
  reactNbrU = reactive({
    NbrU = nrow(reactGDF()%>%filter(urydil == "U"))
  })
  reactNbrnU = reactive({
    NbrnU = nrow(reactGDF()%>%filter(urydil == "/"))
  })
  reactNbrnAdd = reactive({
    NbrnAdd = nrow(reactGDF()) - reactNbrU() - reactNbrnU()
  })
  reactUperc = reactive({
    Uperc = as.numeric(format(round((reactNbrU()/ nrow(reactGDF())) * 100, 2)))
  })
  reactnUperc = reactive({
    nUperc = as.numeric(format(round((reactNbrnU()/ nrow(reactGDF())) * 100, 2)))
  })
  reactnAddperc = reactive({
    nAddperc = as.numeric(format(round(100 - reactUperc() - reactnUperc(), 2)))
  })
  output$plott = renderPlot({
    ggplot(reactGDF()) +
      geom_bar(aes(x = mRNA, fill = reactGDF()$urydil), position = "fill") + 
      scale_y_continuous(breaks = seq(0, 1, 0.1),labels = function(x) paste0(x*100, "%")) +
      geom_text(x = 1, y = as.numeric((reactnAddperc()/2)/100), label = paste(reactnAddperc(), "%"), size = 7*(reactnAddperc()/(reactnAddperc()+1))) +
      geom_text(x = 1, y = as.numeric((reactnAddperc()/100) + ((reactUperc()/2)/100)), label = paste(reactUperc(), "%"), size = 7*(reactUperc()/(reactUperc()+1))) +
      geom_text(x = 1, y = as.numeric((reactnAddperc()/100) + (reactUperc()/100) + ((reactnUperc()/2)/100)), label = paste(reactnUperc(), "%"), size = 7*(reactnUperc()/(reactnUperc()+1))) +
      geom_text(x = 1.15, y = as.numeric((reactnAddperc()/100) + ((reactUperc()/2)/100)), label = paste(reactNbrU(), "read(s)"), size = 5*(reactUperc()/(reactUperc()+1))) +
      geom_text(x = 1.15, y = as.numeric((reactnAddperc()/100) + (reactUperc()/100) + ((reactnUperc()/2)/100)), label = paste(reactNbrnU(), "read(s)"), size = 5*(reactnUperc()/(reactnUperc()+1))) +
      geom_text(x = 1.15, y =as.numeric((reactnAddperc()/2)/100) , label = paste(reactNbrnAdd(), "read(s)"), size = 5*(reactnAddperc()/(reactnAddperc()+1))) +
      scale_fill_discrete(name = "urydilation", labels = c("No uridylation (<60%)", "Urydilation (>60%)", "No additional tail")) +
      theme(legend.title = element_text(size = 18), legend.text = element_text(size = 14)) +
      xlab("Gene") + ylab("distribution in different mRNA")
  })
  observe({
    xrangevall = input$xrangval
    updateSliderInput(session, "xrange", value = xrangevall)
  })
  observe({
      genevall = input$geneval
      updateSelectInput(session, "gene", selected = genevall)
    })
    reactfbindlist <- reactive({
    fbindlist%>%filter(origin==input$barcode)
    })
    output$plot = renderPlot({
    ggplot(reactfbindlist(), aes(x = polya_length, color = chr, )) + geom_density() + xlim(0, input$xrange) +
        scale_color_hue(name = 'Chromosome : ', labels = c("Chromosome 1", "Chromosome 2", "Chromosome 3", "Chromosome 4", "Chromosome 5", "Chloroplastic Chromosome", "Mitochondrial Chromosome")) +
        theme(legend.title = element_text(size = 16), legend.text = element_text(size = 12)) +
        xlab("PolyA tail length") + ylab("Density")
    })
    
    reactTable = reactive({
      intervalTable = rowid_to_column(fbindlist)%>%
        filter(rowid >= input$intervalReads[1], rowid <= input$intervalReads[2])%>%
        separate(col = read_core_id, into = c("id", "chromosome", "mapping_start", "mapping_end"), sep = ",", remove = TRUE)%>%
        mutate(read_length = as.integer(mapping_end) - as.integer(mapping_start),.after = "mapping_end")
    })
    
    output$table = renderDataTable(reactTable(), options = list(scrollX = TRUE))
    
    output$plotreads = renderPlot({
      ggplot(reactTable(), aes(x = read_length, color=type)) + geom_density()
    })
    
    observeEvent(input$intervalReads,{
      newValuesInterval = input$intervalReads
      
      if(valuesInterval[1] != newValuesInterval[1])
        updateSliderInput(session, "intervalReads", value = c(newValuesInterval[1], newValuesInterval[1] + INTERVAL))
      
      if(valuesInterval[2] != newValuesInterval[2])
        updateSliderInput(session, "intervalReads", value = c(newValuesInterval[2] - INTERVAL, newValuesInterval[2]))
      
      valuesInterval <<- newValuesInterval
    })
    
    output$downloadtable = downloadHandler(
      filename = function(){
        paste0("datatable:", valuesInterval[1], "-", valuesInterval[2], ".csv")},
      content = function(file){
        write.csv(reactTable(), file)}
    )
  
}
  
shinyApp(ui, server)
