library(shiny)
library(shinydashboard)
library(tidyverse)
library(here)
library(DT)
library(plotly)
library(leaflet)
library(leaflet.extras)
library(sf)

# Download necessary data -------------------------------------------------
if(!file.exists(here("www", "clover_data.rds"))) {
    viruses <- read_csv(gzcon(url("https://github.com/viralemergence/clover/raw/main/clover/clover_1.0_allpathogens/CLOVER_1.0_Viruses_AssociationsFlatFile.csv")), show_col_types = FALSE)
    
    arenaviridae_hantaviridae <- viruses %>%
        filter(str_detect(PathogenFamily, "arenaviridae|hantaviridae")) %>%
        filter(str_detect(HostOrder, "rodentia")) %>%
        select(AssocID, Host, HostGenus, Pathogen, DetectionMethod, DatabaseVersion, DatabaseDOI, HostTaxID, PathogenTaxID) %>%
        mutate(Reference = paste0("<a href='", DatabaseDOI, "'>", DatabaseVersion, "</a>"),
               Host = str_to_sentence(Host),
               HostGenus = str_to_sentence(HostGenus),
               Pathogen = str_to_sentence(Pathogen)) %>%
        group_by(Host, Pathogen, DetectionMethod) %>%
        arrange(AssocID) %>%
        slice(1) %>%
        select("Association ID" = AssocID,
               "Host species" = Host,
               "Host genus" = HostGenus,
               Pathogen,
               "Detection method" = DetectionMethod,
               "Host ID" = HostTaxID,
               "Pathogen ID" = PathogenTaxID,
               "Database reference" = DatabaseVersion,
               "Database DOI" = DatabaseDOI,
               Reference)
    
    write_rds(arenaviridae_hantaviridae, here("www", "clover_data.rds"))
    
}

clover_data <- read_rds(here("www", "clover_data.rds")) %>%
    ungroup()


# Load and process IUCN data ----------------------------------------------
# Data has been downloaded from the IUCN redlist for all Rodentia, it is then filtered to those species in the clover_data subset
# Names have been manually matched for the 16 that are different in the IUCN database
if(!file.exists(here("www", "IUCN_data", "iucn_ranges.rds"))) {
    unmatched_rodent_names <- tibble(clover_name = unique(clover_data %>%
                                                              filter(!`Host species` %in% IUCN$SCI_NAME) %>%
                                                              pull(`Host species`)),
                                     IUCN_name = c("Necromys amoenus", NA, "Alexandromys fortis",
                                                   "Microtus mystacinus", "Alexandromys maximowiczii",
                                                   "Alexandromys oeconomus", "Clethrionomys gapperi",
                                                   "Clethrionomys glareolus", "Craseomys regulus",
                                                   "Craseomys rufocanus", "Clethrionomys rutilus",
                                                   "Otomys unisulcatus", "Rattus tanezumi", "Neotamias amoenus",
                                                   "Neotamias minimus", "Microtus subterraneus"))
    included_rodents <- clover_data %>%
        left_join(unmatched_rodent_names, by = c("Host species" = "clover_name")) %>%
        mutate(IUCN_name = coalesce(IUCN_name, `Host species`)) %>%
        distinct(`Host species`, IUCN_name)
    
    IUCN <- read_sf(here("www", "IUCN_data", "data_0.shp")) %>%
        filter(SCI_NAME %in% included_rodents$IUCN_name) %>%
        select(SCI_NAME, geometry)
    
    write_rds(IUCN, here("www", "IUCN_data", "iucn_ranges.rds"))
    write_rds(included_rodents, here("www", "IUCN_data", "iucn_match.rds"))
    
} else {
    
    IUCN <- read_rds(here("www", "IUCN_data", "iucn_ranges.rds"))
    included_rodents <- read_rds(here("www", "IUCN_data", "iucn_match.rds"))
}


# UI ----------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- dashboardPage(
    
    skin = "yellow",
    
    dashboardHeader(title = "Arenaviruses and Hantaviruses of rodents",
                    titleWidth = 450),
    
    dashboardSidebar(
        width = 350,
        sidebarMenu(
            menuItem("Homepage", tabName = "homepage", icon = icon("house")),
            menuItem("Known Arenaviridae/Hantaviridae of rodents", tabName = "known_pathogens", icon = icon("viruses")),
            menuItem("Currently known rodent hosts", icon = icon("chart-simple"), startExpanded = TRUE,
                     menuSubItem("Host-Pathogen associations", tabName = "known-h-p", icon = icon("chart-simple")),
                     menuSubItem("Host ranges", tabName = "known-host-ranges", icon = icon("map-location-dot")))
        )
    ),
    
    
    ## Homepage ----------------------------------------------------------------
    
    dashboardBody(
        tabItems(
            # Homepage tab content
            tabItem(tabName = "homepage",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Introduction"),
                            p("Rodents are global hosts of zoonotic pathogens and potential hosts of novel pathogens of epidemic potential. Existing efforts to catalogue host-pathogen associations in these species are limited by global datasets which lack temporal and geographic specificity. Current research is hindered by spatial-, host taxa- and temporal biases within these datasets that are challenging to quantify. Here, we produce a database of studies on rodent-pathogen associations on a global scale focussed on two important viral families, ", strong("Arenaviridae"), " and ", strong("Hantaviridae"), "."),
                            p("Arenaviridae and Hantaviridae are two important, globally distributed, rodent-associated viral families in the order ", strong("Bunyavirales"), ". This order contains several known zoonoses. Arenaviridae associated zoonoses include ", em("Lassa mammarenavirus"), " the cause of Lassa fever in West Africa, ", em("Lymphocytic choriomeningitis virus"), " the cause of Lymphocytic choriomeningitis disease. Hantaviridae associated zoonoses include, ", em("Hantaan orthohantavirus"), " and ", em("Sin Nombre orthohantavirus"), " which cause Hantavirus haemorrhagic fever with renal syndrome and Hantavirus pulmonary syndrome respectively."),
                            p("This data contained within this application has been extracted from a systematic review of the peer-reviewed scientific literature, pre-print articles, ecological reports and `grey` literature with data extracted for subsequent analysis. Data has been extracted on rodent sampling, including species sampled, sampling locations and sampling effort, which is matched to pathogen assays. Where possible these assays have been linked to pathogen sequences archived on NCBI GenBank and pathogenesis studies conducted on these viruses."),
                            p("These data are made available in several formats:"),
                            # Numbered bullet points for the data contained within the application from output$pageslist
                            uiOutput("pageslist"),
                            p(HTML("The repository for the data extraction project and the study protocol is available on GitHub within the <a href=`https://github.com/DidDrog11/arenavirus_hantavirus`>arenavirus_hantavirus</a> project."))
                        )
                    )),
            
            
            ## Known pathogens ---------------------------------------------------------
            
            # Known Arenavirus and Hantavirus pathogens of rodents
            tabItem(tabName = "known_pathogens",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Datasource"),
                            p("CLOVER is a harmonised mammal-virus association database constructed to reconcile and aggregate four datasets that catalogue host-virus associations", HTML("<a href=`https://doi.org/10.1093/biosci/biab080`>(Gibb et al. 2021)</a>"), ". The table below shows all of the contained host-pathogen associations for species of Rodentia and Arenaviridae or Hantaviridae.")),
                        
                        box(width = 12,
                            h1("Table of host-pathogen associations for Rodentia and Arenaviridae or Hantaviridae"),
                            p("151 rodent host species of 48 pathogens are contained within the CLOVER database. 24 of these pathogens are Arenaviridae and 24 Hantaviridae."),
                            downloadButton("download_clover", "Download table as .csv"),
                            br(),
                            p("This will download the entire dataset or a subset depending on the results of the search terms in the search box."),
                            # This datatable contains the output of the clover_data rodent-arenaviridae/hantaviridae subset
                            DTOutput("clover_data"))
                    )),
            
            
            ## Rodent ranges -----------------------------------------------------------
            
            # Plot of rodents known to be hosts of Arenaviridae and Hantaviridae
            tabItem(tabName = "known-h-p",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Known rodent hosts of Arenaviridae and Hantaviridae - CLOVER"),
                            p("Within the CLOVER database 151 rodent species are identified are identified to be hosts of pathogens. Most rodent species are hosts of a single pathogen (84, 56%), with few species host to 5 or more pathogen species (5%)."),
                        ),
                        
                        box(width = 12,
                            height = 700,
                            p("This is an interactive plot of known rodent hosts of Arenaviridae and Hantaviridae. The y-axis shows Host species genera, each bordered section of the bar corresponds to a single species within this genus. Hovering over the bordered section will display the species name and the names of the pathogens associated with that species. The x-axis counts the total number of pathogens within the genus, if multiple host species are hosts of the same pathogen these will be counted multiple times. Not all genera are listed on the y-axis at the default zoom. The plot can be zoomed in on by draggin within the plot area. Double click to resize the plot to its default."),
                            plotlyOutput("rodent_hosts"))
                        
                    )),
            # Ranges of rodents known to be hosts of Arenaviridae and Hantaviridae
            tabItem(tabName = "known-host-ranges",
                    fluidRow(
                        
                        box(width = 4,
                            h1("Filter Host or Pathogen to explore the known range of the rodent host species"),
                            uiOutput("pathogen"),
                            uiOutput("genus"),
                            uiOutput("species")),
                        
                        box(width = 8,
                            h1("Select species to map their rodent range"),
                            DTOutput("table_subset")),
                        
                        box(width = 12,
                            h1("Map of rodent hosts' native range"),
                            p("Rodent range maps have been obtained from the IUCN and NatureServe, the International Union for Conservation of Nature Red List of Threatened Species (2022), downloaded on 2023-05-18. These maps show the native distribution of the selected species in the above table. ", em("n.b. Cavia procellus"), " the Guinea pig does not have a wild distribution and ", em("Rattus flavipectus"), " is classified as a subspecies of ", em("Rattus tanezumi"), ".", strong("Importantly, the pathogens listed on the popups are unlikely to be distributed throughout the range of the rodent host.")),
                            leafletOutput("species_map"))
                        
                    ))
            
        )
    )
    
)

# Server ------------------------------------------------------------------

# Define server logic required to draw a histogram
server <- function(input, output) {
    

# Homepage ----------------------------------------------------------------

    # Numbered bullet points for the data contained within the application
    output$pageslist <- renderUI({
        tags$ol(
            tags$li("The included studies with citations and links to the original publication from which data were extracted."),
            tags$li("A description of the rodent species sampled for pathogens of these two viral families."),
            tags$li("The location of studies and sampling of rodents."),
            tags$li("The prevalence of pathogens of these viral families among sampled rodent species."),
            tags$li("The location of detections of these viruses within hosts."),
            tags$li("Locations of available viral sequences for detected pathogens.")
        )
    })
    

# Known pathogens ---------------------------------------------------------
    
    # Render a datatable of the clover dataframe that is stored in the www folder and loaded in as clover_data
    output$clover_data <- renderDT(
        datatable(clover_data %>%
                      select("Association ID",
                             "Host species",
                             Pathogen,
                             "Detection method",
                             "Host ID",
                             "Pathogen ID",
                             Reference),
                  escape = FALSE)
    )
    
    # This produces the filtered data based on the search field
    output$filtered_clover <- renderPrint({
        input[["clover_data_rows_all"]]
    })
    
    # This button downloads the filtered data
    output$download_clover <- downloadHandler(filename = "clover_rodent_pathogen.csv",
                                              content = function(file) {
                                                  write.csv(clover_data[input[["clover_data_rows_all"]], ] %>%
                                                                select(-Reference),
                                                            file,
                                                            row.names = FALSE)
                                                  })
    

# Rodent ranges -----------------------------------------------------------
    
    output$rodent_hosts <- renderPlotly({
        
        plot_rodent <- clover_data %>%
            ungroup() %>%
            distinct(`Host species`, `Host genus`, Pathogen) %>%
            group_by(`Host species`) %>%
            mutate(`N pathogens` = n()) %>%
            ungroup() %>%
            mutate(`Host genus` = fct_rev(fct_infreq(`Host genus`))) %>%
            group_by(`Host species`) %>%
            mutate(Pathogens = paste0(Pathogen, collapse = ", ")) %>%
            ungroup() %>%
            distinct(`Host species`, `Host genus`, `N pathogens`, Pathogens)
        
        plot_rodent_list <- split(plot_rodent %>%
                                      rename(host_species = `Host species`), seq_len(nrow(plot_rodent)))
        
        plot_ly(x = plot_rodent$`N pathogens`, 
                y = plot_rodent$`Host genus`, 
                stroke = I("black"),
                colors = "YlOrRd",
                customdata = plot_rodent_list, 
                hovertemplate = "Host species: %{customdata.host_species}\nPathogen(s): %{customdata.Pathogens}<extra></extra>",
                type = "bar",
                height = 550) %>%
            layout(barmode = "stack",
                   xaxis = list(title = "Number of pathogens"),
                   yaxis = list(title = "Host genus"))
        
    })
    
    df <- clover_data %>%
        ungroup() %>%
        distinct(`Host genus`, `Host species`, Pathogen)
    
    data <- df
    
    output$table <- renderDataTable({
        if(is.null(data)){return()}
        datatable(data, options = list(scrollX = T))
    })
    
    output$pathogen <- renderUI({
        selectInput(inputId = "Pathogen", "Select pathogen",choices = var_pathogen(), multiple = T)
    })
    output$genus <- renderUI({
        selectInput(inputId = "Genus", "Select genus",choices = var_genus(), multiple = T)
    })
    output$species <- renderUI({
        selectInput(inputId = "Species", "Select species",choices = var_species(), multiple = T)
    })
    
    # Filtered data
    data_filtered <- reactive({
        filter(df,
               Pathogen %in% pathogen(), `Host genus` %in% genus(), `Host species` %in% species())
    })
    
    # Get filters from inputs
    pathogen <- reactive({
        if (is.null(input$Pathogen)) unique(clover_data$Pathogen) else input$Pathogen
    })
    
    genus <- reactive({
        if (is.null(input$Genus)) unique(clover_data$`Host genus`) else input$Genus
    })
    
    species <- reactive({
        if (is.null(input$Species)) unique(clover_data$`Host species`) else input$Species
    })
    
    # Get available categories
    var_pathogen <- reactive({
        file1 <- data
        if(is.null(data)){return()}
        as.list(unique(file1$Pathogen))
    })
    
    var_genus <- reactive({
        filter(data, Pathogen %in% pathogen()) %>% 
            pull(`Host genus`) %>% 
            unique()
    })
    
    var_species <- reactive({
        filter(data, Pathogen %in% pathogen(), `Host genus` %in% genus()) %>% 
            pull(`Host species`) %>% 
            unique()
    })
    
    output$table_subset <- renderDataTable({
        datatable(data_filtered(), options = list(scrollX = T))
    })
    
    output$species_map <- renderLeaflet({
        
        species = data_filtered()[input$table_subset_rows_selected, c("Host species")] %>%
            pull(`Host species`)
        
        iucn_range <- IUCN %>%
            filter(SCI_NAME %in% species) %>%
            rename("Host species" = SCI_NAME)
        
        clover_data_map <- clover_data %>%
            filter(`Host species` %in% species) %>%
            group_by(`Host species`) %>%
            mutate(Pathogens = paste0(Pathogen, collapse = ", ")) %>%
            distinct(`Host species`, Pathogens)
        
        map_species <- left_join(iucn_range, clover_data_map, by = c("Host species"))
        
        leaflet(map_species) %>%
            addPolygons(weight = 1,
                        color = "black",
                        fillColor = "darkred",
                        popup = ~paste0("Rodent species: <em>", `Host species`, "</em><br>Pathogens: ", Pathogens)) %>%
            addProviderTiles("CartoDB.Positron")
        
    })

    
}


# Run application ---------------------------------------------------------

# Run the application 
shinyApp(ui = ui, server = server)
