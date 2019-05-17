###This shiny app creates a table of vegetation suitability and produces a somewhat decent edatopic grid
####Kiri Daust

.libPaths("E:/R packages")
library(shiny)
library(reshape)
library(ggplot2)
library(shinyWidgets)
library(standardize)
library(RColorBrewer)

vegData <- read.csv("BEC10forR.csv")
vegData$SiteUnit <- gsub("[[:space:]]","",vegData$SiteUnit)
vegData$SiteUnit <- gsub('\\$.*',"",vegData$SiteUnit)
vegData$SiteUnit <- gsub('\\..*',"",vegData$SiteUnit)
vegData$SiteUnit <- gsub('/[[:alpha:]]',"/", vegData$SiteUnit)
vegData$SiteUnit <- gsub('/[[:alpha:]]',"/", vegData$SiteUnit)
vegData$SS <- gsub('.*/',"",vegData$SiteUnit)
vegData$BGC <- gsub("/.*","",vegData$SiteUnit)

vegData[is.na(vegData)] <- 0
vegData$Cover <- vegData$TotalA + vegData$TotalB +vegData$Cover6 + vegData$Cover7
vegSpec <- as.character(unique(vegData$Species[vegData$Lifeform == 3 | vegData$Lifeform == 4]))
vegSpec <- vegSpec[order(vegSpec)]
vegData <- vegData[,c(2,10,9,3,11)]

edatope <- read.csv("Edatope_v10.csv")
edatope$numeric <- edatope$Edatopic
edatope$Edatopic <- gsub("[[:alpha:]]","", edatope$Edatopic)
edatope$numeric <- gsub("[[:digit:]]","", edatope$numeric)
edatope <- edatope[,c(2,3,6,4)]
colnames(edatope) <- c("BGC","Unit","Alpha","Numeric")

ui <- fluidPage(
   
   # Application title
   titlePanel("Veg Suitability"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
          htmlOutput("selectBGC"),
          pickerInput("Species",
                      "Choose Species:",
                      choices = vegSpec,
                      multiple = TRUE
        )
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         tableOutput("Table"),
         plotOutput("Plot")
      )
   )
)

# Define server logic 
server <- function(input, output) {
  
  getSuit <- reactive({
    notSuit <- as.character(unique(vegData$SiteUnit))
    vaccDat <- vegData[vegData$Species == input$Species[1],]
    vaccDat <- vaccDat[order(vaccDat$SiteUnit),]
    vaccDat$SiteUnit <- as.character(vaccDat$SiteUnit)
    vaccDat$SiteUnit <- as.factor(vaccDat$SiteUnit)
    vaccDat$Mean <- ave(vaccDat$Cover, vaccDat$SiteUnit, FUN = mean)
    vaccDat <- unique(vaccDat[,c(4,1,6)])
    vaccDat$SS <- gsub('.*/',"",vaccDat$SiteUnit)
    vaccDat$BGC <- gsub("/.*","",vaccDat$SiteUnit)
    vaccDat$Std <- scale(vaccDat$Mean)
    vaccDat <- vaccDat[complete.cases(vaccDat),]
    vaccDat$Suit <- ifelse(vaccDat$Std <= -0.8, 3,
                           ifelse(vaccDat$Std > -0.8 & vaccDat$Std <= 0.8, 2,
                                  ifelse(vaccDat$Std > 0.8, 1, NA)))
    
    
    notSuit <- notSuit[!(notSuit %in% vaccDat$SiteUnit)]
    notSuit <- data.frame(SiteUnit = notSuit, Suit = rep(0,length(notSuit)), 
                          Species = rep(input$Species[1],length(notSuit)), 
                          Mean = rep(0,length(notSuit)), Std = rep(0,length(notSuit)))
    notSuit$SS <- gsub('.*/',"",notSuit$SiteUnit)
    notSuit$BGC <- gsub("/.*","",notSuit$SiteUnit)
    vaccDat <- rbind(vaccDat, notSuit)
    colnames(vaccDat)[2] <- "Unit"
    
    if(length(input$Species) > 1){
      for (i in 2:length(input$Species)){
        notSuit <- as.character(unique(vegData$SiteUnit))
        temp <- vegData[vegData$Species %in% input$Species[i],]
        temp <- temp[order(temp$SiteUnit),]
        temp$SiteUnit <- as.character(temp$SiteUnit)
        temp$SiteUnit <- as.factor(temp$SiteUnit)
        temp$Mean <- ave(temp$Cover, temp$SiteUnit, FUN = mean)
        temp <- unique(temp[,c(4,1,6)])
        temp$SS <- gsub('.*/',"",temp$SiteUnit)
        temp$BGC <- gsub("/.*","",temp$SiteUnit)
        temp$Std <- scale(temp$Mean)
        temp <- temp[complete.cases(temp),]
        temp$Suit <- ifelse(temp$Std <= -0.8, 3,
                               ifelse(temp$Std > -0.8 & temp$Std <= 0.8, 2,
                                      ifelse(temp$Std > 0.8, 1, NA)))
        
        
        notSuit <- notSuit[!(notSuit %in% temp$SiteUnit)]
        notSuit <- data.frame(SiteUnit = notSuit, Suit = rep(0,length(notSuit)), 
                              Species = rep(input$Species[i],length(notSuit)), 
                              Mean = rep(0,length(notSuit)), Std = rep(0,length(notSuit)))
        notSuit$SS <- gsub('.*/',"",notSuit$SiteUnit)
        notSuit$BGC <- gsub("/.*","",notSuit$SiteUnit)
        temp <- rbind(temp, notSuit)
        colnames(temp)[2] <- "Unit"
        vaccDat <- rbind(vaccDat,temp)
      }
    }
    return(vaccDat)
  })
   
   output$Table <- renderTable({
     if(length(input$Species) >= 1){
     vaccDat <- getSuit()
     vegSub <- vaccDat[vaccDat$BGC == input$BGC,]
     tableOut <- cast(vegSub, Species ~ SS, value = "Suit", fun.aggregate = mean)
     }
   })
   
   output$selectBGC <- renderUI({
     vaccDat <- getSuit()
     pickerInput("BGC","Choose Zone", choices = as.character(unique(vaccDat$BGC[vaccDat$Suit != 0])))
   })
   
   output$Plot <- renderPlot({
     if(length(input$Species == 1)){
     vaccDat <- getSuit()
     vaccDat <- vaccDat[,-1]
     vegSub <- vaccDat[vaccDat$BGC == input$BGC,]
     edaSub <- edatope[edatope$BGC == input$BGC,]
     edaSub <- merge(vegSub, edaSub, by = "Unit", all = TRUE)
     edaSub <- edaSub[complete.cases(edaSub),c(8,9,6,3)]
     edaSub$SS <- as.factor(edaSub$SS)
     edaSub$Numeric <- as.numeric(edaSub$Numeric)
     edaSub$Numeric <- edaSub$Numeric * -1
     edaSub$Suit <- as.numeric(edaSub$Suit)
     edaSub$xmin <- ave(edaSub$Alpha, edaSub$SS, FUN = min)
     edaSub$xmax <- ave(edaSub$Alpha, edaSub$SS, FUN = max)
     edaSub$ymin <- ave(edaSub$Numeric, edaSub$SS, FUN = min)
     edaSub$ymax <- ave(edaSub$Numeric, edaSub$SS, FUN = max)
     
     edaGrid <- edaSub[,4:8]
     edaGrid <- unique(edaGrid)
     edaGrid$ymax <- as.numeric(edaGrid$ymax)
     edaGrid$xmint <- lapply(edaGrid$xmin, function(x){which(LETTERS == x)}) 
     edaGrid$xmaxt <- lapply(edaGrid$xmax, function(x){which(LETTERS == x)}) 
     edaGrid$xmint <- as.numeric(edaGrid$xmint)
     edaGrid$xmaxt <- as.numeric(edaGrid$xmaxt)
     edaGrid$ymin <- as.numeric(edaGrid$ymin)
     edaGrid$xmint[edaGrid$xmint != 1] <- jitter(edaGrid$xmint[edaGrid$xmint != 1], amount = 0.1)
     edaGrid$xmaxt[edaGrid$xmaxt != max(edaGrid$xmaxt)] <- jitter(edaGrid$xmaxt[edaGrid$xmaxt != max(edaGrid$xmaxt)], amount = 0.1)
     edaGrid$ymin[edaGrid$ymin != min(edaGrid$ymin)] <- jitter(edaGrid$ymin[edaGrid$ymin != min(edaGrid$ymin)], amount = 0.1)
     edaGrid$ymax[edaGrid$ymax != -1] <- jitter(edaGrid$ymax[edaGrid$ymax != -1], amount = 0.1)
     edaSub$Numeric <- as.factor(edaSub$Numeric)
     edaSub$Alpha <- lapply(edaSub$Alpha, function(x){which(LETTERS == x)})
     edaSub$Alpha <- as.numeric(edaSub$Alpha)
     edaSub$Numeric <- as.numeric(as.character(edaSub$Numeric))
     
     
     edaSub$Suit <- as.factor(edaSub$Suit)
     edaSub <- edaSub[edaSub$Suit != 0,]
     ggplot(data = edaSub) +
       geom_raster(aes(x= Alpha, y = Numeric, fill = Suit))+
       scale_fill_manual(values = c("green", "orange","red"))+
       geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5), color = "grey", size = 0.8)+
       geom_hline(yintercept = seq(-0.5,-7.5,-1), color = "grey", size = 0.8)+
       geom_rect(data = edaGrid, fill = "NA", col = rainbow(length(edaGrid$SS)), aes(xmin = xmint-0.5, xmax = xmaxt+0.5, ymin = ymin-0.5, ymax = ymax+0.5), size = 1.5)+
       annotate("text", x = edaGrid$xmint, y = edaGrid$ymax, label = edaGrid$SS, col = "black") +
       scale_x_continuous(breaks = c(1,2,3,4,5), labels = c("A","B","C","D","E"), position = "top")+
       scale_y_continuous(breaks = c(-1:-8), labels = c(1,2,3,4,5,6,7,8))+
       theme_classic()+
       coord_fixed()+ xlab("Soil Nutrient Regime (poor to rich)") +
       ylab("Soil Moisture Regime")
   }
       
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

