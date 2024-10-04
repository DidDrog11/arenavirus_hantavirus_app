library(shiny)
library(shinydashboard)
library(shinyjs)
library(tidyverse)
library(here)
library(DT)
library(plotly)
library(leaflet)
library(leaflet.extras)
library(sf)
library(cowplot)
library(bib2df)
library(bibtex)

# Download CLOVER data -------------------------------------------------
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


# Download review data ----------------------------------------------------
update_data = FALSE

if(update_data) {
    review_data <- readRDS(gzcon(url("https://github.com/DidDrog11/arenavirus_hantavirus/raw/main/data/clean_data/2024-10-04_data.rds")))
    
    names(review_data$pathogen)[names(review_data$pathogen) == "virus_clean"] <- "scientificName"
    names(review_data$pathogen)[names(review_data$pathogen) == "n_positive"] <- "occurrenceRemarks"
    
    write_rds(review_data, here("www", "review_data.rds"))
    
} else {
    
    review_data <- read_rds(here("www", "review_data.rds"))
}

# Format and rename citation columns --------------------------------------

citations <- review_data$citations %>%
    drop_na(study_id) %>%
    select(study_id, AUTHOR = Author, YEAR = `Publication Year`, TITLE = Title, JOURNAL = `Publication Title`, ISSN, ISSUE = Issue, VOLUME = Volume, DOI) %>%
    mutate(AUTHOR = lapply(strsplit(AUTHOR, ";"), function(x) trimws(x)),
           ISSUE = as.numeric(ISSUE),
           VOLUME = as.numeric(VOLUME),
           CATEGORY = case_when(!is.na(JOURNAL) ~ "Article",
                                TRUE ~ NA)) %>%
    rowwise() %>%
    mutate(BIBTEXKEY = paste0(str_remove_all(str_split(AUTHOR[[1]][1], ",", simplify = TRUE)[, 1], " "), YEAR, "_", study_id))


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
                     menuSubItem("Host ranges", tabName = "known-host-ranges", icon = icon("map-location-dot"))),
            menuItem("Included studies", tabName = "studies", icon = icon("book")),
            menuItem("Rodent sampling locations", icon = icon("chart-simple"), startExpanded = TRUE,
                     menuSubItem("Sampling locations", tabName = "sampling-locations", icon = icon("chart-simple")),
                     menuSubItem("Sampling locations by rodent/pathogen species", tabName = "sampling-locations-stratified", icon = icon("map-location-dot")),
                     menuSubItem("Sequencing locations by rodent/pathogen species", tabName = "sequence-locations", icon = icon("dna")))
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
                            p(em("n.b. Throughout, where the term rodent is used I include all small mammals. This is incorrect and will be addressed as the project nears completion. The most common taxa that will be grouped in this way are shrews which are not rodents.")),
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
                            p("Rodent range maps have been obtained from the IUCN and NatureServe, the International Union for Conservation of Nature Red List of Threatened Species (2022), downloaded on 2023-05-18. These maps show the native distribution of the selected species in the above table. ", em("n.b. Cavia procellus"), " the Guinea pig does not have a wild distribution and ", em("Rattus flavipectus"), " is classified as a subspecies of ", em("Rattus tanezumi"), ".", strong("Importantly, the pathogens listed on the popups are unlikely to be distributed throughout the entire range of the rodent host.")),
                            leafletOutput("species_map"))
                        
                    )),
            

            ## Included studies ----
            tabItem(tabName = "studies",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Included studies"),
                            p("Following a systematic search of the available literature and review for suitability of inclusion, we have identified ", strong(paste(nrow(review_data$citations %>% filter(!str_detect(decision, "Exclude"))))),
                              " studies for review of full texts and data extraction."),
                            p("As of 2024-01-30 data have been extracted from ", 
                              strong(paste(nrow(review_data$citations %>% filter(!is.na(study_id))))),
                              paste0("(", round(nrow(review_data$citations %>% filter(!is.na(study_id)))/nrow(review_data$citations %>% filter(str_detect(decision, "Include"))) * 100, 0), "%)"),
                              " studies.",  "This equates to ", paste0(round(nrow(review_data$citations %>% filter(!is.na(study_id)))/nrow(review_data$citations %>% filter(str_detect(processed, "y"))) * 100, 0), "%"), "of studies selected for full-text review (that have been reviewed), having data that meets inclusion and exclusion criteria."),
                            plotOutput("review_status")),
                        
                        box(width = 12,
                            h1("Details of included studies"),
                            DTOutput("included_studies")),
                        
                        box(width = 12,
                            h1("Timeline of included studies"),
                            p("Effort to investigate the prevalence of Arenaviruses and Hantaviruses has changed over time, this can be seen through both the number of entries within NCBI PUBMED for the search terms `rodent*` AND (`arenavir*` OR `hantavir*`) and the publication dates of included studies."),
                            plotOutput("included_studies_timeline"))
                    )),
            
            ## Sampled rodents ----
            tabItem(tabName = "sampling-locations",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Sampling locations"),
                            p("Data has been extracted on the locations of rodent species that have been detected (typically through rodent trapping). The date during which sampling occurred is reported at the highest resolution obtainable from the study. The point relates to the position at which the rodent was detected, either the trap, trapping grid, study area or region. The locations are dependent on the level of detail given in the original study. Selecting a point will provide additional information including; whether the species was present or absent at that location, the number of individuals detected, the locality, the habitat type in which the trap was placed and coordinate resolution of the point."),
                            p("The map below shows the locations of sampling contained within this dataset. Selecting sampling points on the map will populate the table with the either the studies containing these sampling locations or the species and pathogens identified at these locations. Both datasets can then be downloaded, formatted as they are in the table."))),
                    
                    fluidRow(
                        
                        box(width = 12,
                            h1("Mapping sampling"),
                            p(""),
                        leafletOutput("sampling_locations",
                                      height = 600)))),
            
            tabItem(tabName = "sampling-locations-stratified",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Sampled rodents stratified by pathogen or rodent species"),
                            p("This map and the associated tables contain the same information as the previous page. Here, data may be filtered based on the rodent species or pathogen species with the filtered subset mapped.",
                              HTML(paste("<ul>
                                         <li>The dataset contains", nrow(review_data$host), "rodent records.</li>
                                         <li>These include data on", length(unique(review_data$host$species)), "rodent species.</li>
                                         <li>There are data on", length(unique(review_data$pathogen$scientificName)), "pathogen species.</li>
                                         <li>There are data on", nrow(review_data$pathogen %>%
                                             filter(occurrenceRemarks >= 1) %>%
                                             distinct(host_species, scientificName)), "host-pathogen associations (including negative associations).</li>
                                         </ul>")))),
                        
                        box(width = 12,
                            h1("Filter data on host or pathogen"),
                            column(3, uiOutput("pathogen_family_rev")),
                            column(3, uiOutput("pathogen_species_rev")),
                            column(3, uiOutput("host_genus_rev")),
                            column(3, uiOutput("host_species_rev")),
                            collapsible = TRUE),
                        
                        box(width = 12,
                            p("Hiding the side bar (three lines next to the title) may help fitting the table. The download button will download the selected data (filtered_data.csv) and associated citations (filtered_citations.bib).\n", 
                              downloadButton("download_path_data", "Download Filtered Data and Citations")),
                            DTOutput("review_table_subset"),
                            collapsible = TRUE),
                        
                        box(width = 12,
                            leafletOutput("pathogen_map", height = 600))
                    )),
            tabItem(tabName = "sequence-locations",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Sequences deposited in NCBI GenBank"),
                            p("This map and the associated tables contain information on sequences available in NCBI GenBank. Metadata is often incomplete to accompany the sequences. Our current data synthesis approach enriches sequence metadata by linking obtained sequences back to sampled rodents. This improves the geographic resolution of sampling locations and dates sampled which may be helpful for interpreting analyses using these sequences.")),
                        
                        box(width = 12,
                            h1("Filter data on host or pathogen"),
                            column(3, uiOutput("seq_pathogen_family_rev")),
                            column(3, uiOutput("seq_pathogen_species_rev")),
                            column(3, uiOutput("seq_host_genus_rev")),
                            column(3, uiOutput("seq_host_species_rev")),
                            collapsible = TRUE),
                        
                        box(width = 12,
                            p("Hiding the side bar (three lines next to the title) may help fitting the table. The download button will download the selected data (filtered_data.csv) and associated citations (filtered_citations.bib).\n", 
                              downloadButton("download_seq_data", "Download Filtered Sequence Data and Citations")),
                            DTOutput("seq_review_table_subset"),
                            collapsible = TRUE),
                        
                        box(width = 12,
                           leafletOutput("sequence_map", height = 600))
                        
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
                   yaxis = list(title = "Host genus")) %>% 
            layout(yaxis = list(tickmode = "auto", nticks = length(unique(plot_rodent$`Host genus`))))
        
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
        selectInput(inputId = "Genus", "Select host genus",choices = var_genus(), multiple = T)
    })
    output$species <- renderUI({
        selectInput(inputId = "Species", "Select host species",choices = var_species(), multiple = T)
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
            mutate(Pathogens = paste0(unique(Pathogen), collapse = ", ")) %>%
            distinct(`Host species`, Pathogens)
        
        map_species <- left_join(iucn_range, clover_data_map, by = c("Host species"))
        
        leaflet(map_species) %>%
            addPolygons(weight = 1,
                        color = "black",
                        fillColor = "darkred",
                        popup = ~paste0("Rodent species: <em>", `Host species`, "</em><br>Pathogens: ", Pathogens)) %>%
            addProviderTiles("CartoDB.Positron")
        
    })
    

# Included studies --------------------------------------------------------

    output$review_status <- renderPlot({
        
        all_studies <- review_data$citations %>%
            nrow()
        
        excluded_full_text <- review_data$citations %>%
            filter(str_detect(decision, "Exclude")) %>%
            nrow()
        
        extracted_full_text <- review_data$citations %>%
            filter(!is.na(study_id)) %>%
            nrow() 
        
        awaiting_review <- all_studies - excluded_full_text - extracted_full_text
        
        tibble(status = c("Extracted", "Awaiting review", "Excluded full text"),
               n_studies = c(extracted_full_text, awaiting_review, excluded_full_text)) %>%
            mutate(status = factor(status, levels = c( "Excluded full text", "Awaiting review", "Extracted")),
                   percentage = round(n_studies / sum(n_studies) * 100, 0)) %>%
            ggplot(aes(x = "", y = n_studies, fill = status)) +
            geom_col(width = 0.7) +  # Use geom_col for a single bar
            geom_text(aes(label = paste0(percentage, "%")), 
                      position = position_stack(vjust = 0.5), 
                      size = 8, color = "white") +
            labs(title = "Data extraction progress",
                 x = NULL,
                 y = "Number of Studies",
                 fill = element_blank()) +
            theme_minimal() +
            theme(axis.text.y = element_blank(),  # Hide x-axis labels
                  legend.position = "top") +
            coord_flip()
        
        })
    
    output$included_studies <- renderDataTable({
        
        
        doi_links <- review_data$citations %>%
            mutate(DOI = paste0("<a href='https://doi.org/", gsub("#", "%23", DOI), "'>", DOI, "</a>")) %>%
            select(full_text_id, DOI, journal = `Publication Title`)
        
        review_data$citations %>%
            filter(!is.na(study_id)) %>%
            select(study_id, full_text_id, Author, Title, `Publication Year`) %>%
            mutate(Author = str_split(Author, ";", simplify = TRUE)[, 1]) %>%
            left_join(doi_links, by = c("full_text_id")) %>%
            select("Study ID" = study_id, 
                   "First author" = Author, "Year" = `Publication Year`,
                   Title, "Journal" = journal, DOI) %>%
            datatable(escape = FALSE, rownames = FALSE)
        
        
    })
    
    output$included_studies_timeline <- renderPlot({
        
        search_results_plot <- review_data$citations %>%
            mutate(year = as.numeric(`Publication Year`)) %>%
            group_by(year) %>%
            summarise(n = n()) %>%
            ggplot() +
            geom_col(aes(x = year, y = n, fill = year)) +
            scale_fill_viridis_c("magma") +
            labs(x = "Year",
                 y = "Number of search results") +
            guides(fill = "none") +
            theme_bw()
        
        included_studies_plot <- review_data$citations %>%
            drop_na(study_id) %>%
            mutate(publication_year = as.numeric(`Publication Year`)) %>%
            drop_na(publication_year) %>%
            group_by(publication_year) %>%
            summarise(n = n()) %>%
            ggplot() +
            geom_col(aes(x = publication_year, y = n, fill = publication_year)) +
            scale_fill_viridis_c("magma", limits = c(min(as.numeric(citations$YEAR), na.rm = TRUE),
                                                     max(as.numeric(citations$YEAR), na.rm = TRUE))) +
            labs(x = "Year",
                 y = "Number of included studies") +
            guides(fill = "none") +
            theme_bw() +
            coord_cartesian(xlim = c(min(as.numeric(citations$YEAR), na.rm = TRUE),
                                     max(as.numeric(citations$YEAR), na.rm = TRUE)))
        
        plot_grid(plotlist = list(search_results_plot,
                                  included_studies_plot),
                  nrow = 2)
        
    })
    

# Rodent locations --------------------------------------------------------

    output$sampling_locations <- renderLeaflet({
        
        locations <- review_data$host %>%
            group_by(study_id, eventDate, locality, country, verbatimLocality, coordinate_resolution, 
                     decimalLatitude, decimalLongitude) %>%
            summarise(p_species = sum(as.numeric(individualCount) >= 1, na.rm = TRUE),
                      a_species = sum(as.numeric(individualCount) == 0), na.rm = TRUE) %>%
            mutate(low_resolution = case_when(is.na(locality) ~ TRUE,
                                              coordinate_resolution == "country" ~ TRUE,
                                              coordinate_resolution == "state" ~ TRUE,
                                              coordinate_resolution == "region" ~ TRUE,
                                              coordinate_resolution == "province" ~ TRUE,
                                              coordinate_resolution == "district" ~ TRUE,
                                              coordinate_resolution == "county region" ~ TRUE,
                                              coordinate_resolution == "county" ~ TRUE,
                                              coordinate_resolution == "municipality" ~ TRUE,
                                              TRUE ~ FALSE)) %>%
            drop_na(decimalLatitude) %>%
            drop_na(decimalLongitude) %>%
            left_join(review_data$citations %>%
                          drop_na(study_id) %>%
                          mutate(Author = str_split(Author, ";", simplify = TRUE)[, 1]) %>%
                          select(study_id, full_text_id, identifiedBy = Author, datasetName = Title, `Publication Year`, DOI), by = "study_id") %>%
            st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
        
        res_color_palette <- c("darkred", "darkgreen")
        res_color_factor <- colorFactor(palette = res_color_palette, domain = locations$low_resolution)
        
        
        leaflet(locations) %>%
            addTiles() %>%
            addCircleMarkers(color = ~res_color_factor(low_resolution),
                             stroke = FALSE,
                             fillOpacity = 0.8,
                             clusterOptions = NULL,
                             radius = 3,
                             popup = ~paste0(
                                 "• Number of rodent species detected: ", p_species, "<br/>",
                                 "• Number of rodent species not detected: ", a_species, "<br/>",
                                 "• Study ID: ", study_id, "<br/>",
                                 "• ", locality, ", ", country, "<br/>",
                                 "• Coordinate resolution: ", str_to_sentence(coordinate_resolution), "<br/>",
                                 "• ", identifiedBy, ", et al. ", `Publication Year`, ". ", datasetName, "<br/>",
                                 "• DOI: ", DOI)
                             ) %>%
            addLegend(title = "Low Coordinate Resolution",
                      colors = res_color_palette,
                      values = unique(locations$low_resolution),
                      labels = c("FALSE", "TRUE"),
                      opacity = 1
                      )
        
    })


    # Stratified records ------------------------------------------------------

    data_rev <- reactiveVal(review_data$pathogen)
    
    output$pathogen_family_rev <- renderUI({
        selectizeInput(
            inputId = "pathogen_family_rev",
            label = "Select pathogen family",
            choices = var_pathogen_family_rev(),
            multiple = TRUE,
            selected = c("Hantaviridae", "Mammarenaviridae"),
            options = list(
                server = TRUE,
                placeholder = 'Select pathogen family'
            )
        )
    })
    output$pathogen_species_rev <- renderUI({
        selectizeInput(
            inputId = "pathogen_species_rev",
            "Select pathogen species",
            choices = var_pathogen_species_rev(),
            multiple = TRUE,
            options = list(
                server = TRUE,
                placeholder = 'Select pathogen species'
            ))
    })
    output$host_genus_rev <- renderUI({
        selectizeInput(
            inputId = "host_genus_rev",
            "Select host genus",
            choices = var_host_genus_rev(),
            multiple = TRUE,
            options = list(
                server = TRUE,
                placeholder = 'Select host genus'
            ))
    })
    output$host_species_rev <- renderUI({
        selectizeInput(
            inputId = "host_species_rev",
            "Select host species",
            choices = var_host_species_rev(),
            multiple = TRUE,
            options = list(
                server = TRUE,
                placeholder = 'Select host species'
            ))
    })
    
    # Get available categories
    var_pathogen_family_rev <- reactive({
        sort(unique(review_data$pathogen$family))
    })
    
    var_pathogen_species_rev <- reactive({
        
        if (!is.null(input$pathogen_family_rev)) {
            sort(data_rev() %>%
                filter(family %in% input$pathogen_family_rev) %>%
                distinct(scientificName) %>%
                pull())
        } else {
            sort(data_rev() %>%
                distinct(scientificName) %>%
                pull())
        }

    })
    
    var_host_genus_rev <- reactive({
        sort(unique(review_data$pathogen$host_genus))
    })
    
    var_host_species_rev <- reactive({
        
        if (!is.null(input$host_genus_rev)) {
            sort(data_rev() %>%
                filter(host_genus %in% input$host_genus_rev) %>%
                distinct(host) %>%
                pull())
        } else {
            sort(data_rev() %>%
                distinct(host) %>%
                pull())
        }
        
    })
    
    # Filtered data
    data_filtered_rev <- reactive({
        filtered_data <- data_rev()
        
        # Filter based on pathogen family
        if (!is.null(input$pathogen_family_rev) & is.null(input$pathogen_species_rev)) {
            filtered_data <- filter(filtered_data, family %in% input$pathogen_family_rev)
        }
        
        # Filter based on pathogen species
        if (!is.null(input$pathogen_species_rev)) {
            filtered_data <- filter(filtered_data, scientificName %in% input$pathogen_species_rev)
        }
        
        # Filter based on host genus
        if (!is.null(input$host_genus_rev) & is.null(input$host_species_rev)) {
            filtered_data <- filter(filtered_data, host_genus %in% input$host_genus_rev)
        }
        
        # Filter based on host species
        if (!is.null(input$host_species_rev)) {
            filtered_data <- filter(filtered_data, host %in% input$host_species_rev)
        }
        
        return(filtered_data)
    })
    
    output$review_table_subset <- renderDataTable({
        datatable(data_filtered_rev(), options = list(scrollX = TRUE, autoWidth = TRUE),
                  # colnames = c("Path ID", "Rodent ID", "Study ID", "Host genus", "Host species", "Locality", "Country", "Coord Resolution",
                  #              "Lat", "Long", "Pathogen rank", "Path family", "Path species", "Assay", "Tested", "Negative", "Positive",
                  #              "Inconclusive", "Note"),
                  class = 'custom-datatable',
                  rownames = FALSE)
    })
    
    # Function to filter citations based on selected study IDs
    filterCitations <- function(selected_study_ids) {
        citations %>%
            filter(study_id %in% selected_study_ids)
    }
    
    # Function to create BibTeX file from filtered citations
    createBibTeXFile <- function(filtered_citations, file_path) {
        if (nrow(filtered_citations) > 0) {
            # Print the filtered_citations for debugging
            print(filtered_citations)
            
            # Convert to BibTeX format
            bib_data <- df2bib(filtered_citations)
            
            # Create BibTeX entries
            bib_entries <- apply(filtered_citations, 1, function(row) {
                entry <- paste0("@Article{", row$BIBTEXKEY, ",\n",
                                "  study_id = {", row$study_id, "},\n",
                                "  Author = {", paste(row$AUTHOR, collapse = "; "), "},\n",
                                "  Year = {", row$YEAR, "},\n",
                                "  Title = {", row$TITLE, "},\n",
                                "  Journal = {", row$JOURNAL, "},\n",
                                "  Issn = {", row$ISSN, "},\n",
                                "  Issue = {", row$ISSUE, "},\n",
                                "  Volume = {", row$VOLUME, "},\n",
                                "  Doi = {", row$DOI, "}\n",
                                "}\n")
                return(entry)
            })
            
            # Write BibTeX entries to the file
            writeLines(bib_entries, con = file_path)
        }
    }
    
    output$download_path_data <- downloadHandler(
        filename = function() {
            paste("filtered_data_and_citations_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file) {
            # Create a temporary directory to store files
            temp_dir <- tempdir()
            
            # Filtered table
            write.csv(data_filtered_rev(), file.path(temp_dir, "filtered_data.csv"))
            
            # Filtered citations
            selected_study_ids <- unique(data_filtered_rev()$study_id)
            filtered_citations <- filterCitations(selected_study_ids)
            bib_path <- file.path(temp_dir, "filtered_citations.bib")
            createBibTeXFile(filtered_citations, bib_path)
            
            # Create a zip file with both CSV and BibTeX files
            zip(file, files = c(file.path(temp_dir, "filtered_data.csv"), bib_path), flags = "-j")
        }
    )

# Plot pathogen records ---------------------------------------------------

    # Initialize an empty map
    output$pathogen_map <- renderLeaflet({
        leaflet() %>%
            addTiles()
    })
    
    colour_palette <- c("darkgreen", "darkred")
    
    # ReactiveVal to track filter status
    filter_applied <- reactiveVal(FALSE)
    
    # Update the map based on filtered data
    observe({
        filtered_data <- data_filtered_rev() %>%
            filter(occurrenceRemarks >= 1)
        
        # Check if there are filters applied
        if ((!is.null(input$pathogen_family_rev) && length(input$pathogen_family_rev) > 0) ||
            (!is.null(input$pathogen_species_rev) && length(input$pathogen_species_rev) > 0) ||
            (!is.null(input$host_genus_rev) && length(input$host_genus_rev) > 0) ||
            (!is.null(input$host_species_rev) && length(input$host_species_rev) > 0)) {
            
            # Filter the data based on selected filters
            filtered_data <- filtered_data %>%
                filter(
                    if (!is.null(input$pathogen_family_rev)) family %in% input$pathogen_family_rev else TRUE,
                    if (!is.null(input$pathogen_species_rev)) scientificName %in% input$pathogen_species_rev else TRUE,
                    if (!is.null(input$host_genus_rev)) host_genus %in% input$host_genus_rev else TRUE,
                    if (!is.null(input$host_species_rev)) associatedTaxa %in% input$host_species_rev else TRUE
                )
            
            # Check if any records are left after filtering
            if (nrow(filtered_data) > 0) {
                # Map organismQuantity to colour
                colours <- ifelse(filtered_data$n_positive >= 1, "darkred", "darkgreen")
                
                # Clear existing markers and add new ones
                leafletProxy("pathogen_map") %>%
                    clearMarkers() %>%
                    addCircleMarkers(data = filtered_data,
                                     lng = ~decimalLongitude,
                                     lat = ~decimalLatitude,
                                     radius = 8,
                                     color = colours,
                                     fillOpacity = 0.7,
                                     popup = ~paste("Study ID: ", study_id, "<br>",
                                                    "Rodent ID: ", associated_rodent_record_id, "<br>",
                                                    "Rodent: ", host, "<br>",
                                                    "N rodents tested: ", n_tested, "<br>",
                                                    "Pathogen ID: ", pathogen_record_id, "<br>",
                                                    "Pathogen: ", scientificName, "<br>",
                                                    "N pathogen detected: ", n_positive),
                                     clusterOptions = markerClusterOptions(
                                         spiderfyDistanceMultiplier = 2
                                     )) %>%
                    fitBounds(
                        lat1 = min(filtered_data$decimalLatitude, na.rm = TRUE),
                        lng1 = min(filtered_data$decimalLongitude, na.rm = TRUE),
                        lat2 = max(filtered_data$decimalLatitude, na.rm = TRUE),
                        lng2 = max(filtered_data$decimalLongitude, na.rm = TRUE)
                    )
                
                # Update filter status
                filter_applied(TRUE)
            } else {
                # If no records match the filters, clear markers
                leafletProxy("pathogen_map") %>%
                    clearMarkers()
                
                # Update filter status
                filter_applied(FALSE)
            }
        } else {
            # If no filters applied, clear markers
            leafletProxy("pathogen_map") %>%
                clearMarkers()
            
            # Update filter status
            filter_applied(FALSE)
        }
    })

    
    # Clear markers when filter is removed
    observe({
        if (!filter_applied()) {
            leafletProxy("pathogen_map") %>%
                clearMarkers()
        }
    })
    


# Stratify sequence records -----------------------------------------------

    seq_data_rev <- reactiveVal(review_data$sequences %>%
                                    mutate(accession_number = sprintf('<a href="https://www.ncbi.nlm.nih.gov/search/all/?term=%s" target="_blank">%s</a>',
                                                                      accession_number,
                                                                      accession_number)))
    
    output$seq_pathogen_family_rev <- renderUI({
        selectizeInput(
            inputId = "seq_pathogen_family_rev",
            label = "Select pathogen family",
            choices = seq_var_pathogen_family_rev(),
            multiple = TRUE,
            selected = c("Hantaviridae", "Arenaviridae", ""),
            options = list(
                server = TRUE,
                placeholder = 'Select pathogen family'
            )
        )
    })
    output$seq_pathogen_species_rev <- renderUI({
        selectizeInput(
            inputId = "seq_pathogen_species_rev",
            "Select pathogen species",
            choices = seq_var_pathogen_species_rev(),
            multiple = TRUE,
            options = list(
                server = TRUE,
                placeholder = 'Select pathogen species'
            ))
    })
    output$seq_host_genus_rev <- renderUI({
        selectizeInput(
            inputId = "seq_host_genus_rev",
            "Select host genus",
            choices = seq_var_host_genus_rev(),
            multiple = TRUE,
            options = list(
                server = TRUE,
                placeholder = 'Select host genus'
            ))
    })
    output$seq_host_species_rev <- renderUI({
        selectizeInput(
            inputId = "seq_host_species_rev",
            "Select host species",
            choices = seq_var_host_species_rev(),
            multiple = TRUE,
            options = list(
                server = TRUE,
                placeholder = 'Select host species'
            ))
    })
    
    # Get available categories
    seq_var_pathogen_family_rev <- reactive({
        sort(unique(review_data$sequences$family))
    })
    
    seq_var_pathogen_species_rev <- reactive({
        
        if (!is.null(input$seq_pathogen_family_rev)) {
            sort(seq_data_rev() %>%
                     filter(family %in% input$seq_pathogen_family_rev) %>%
                     distinct(virus_clean) %>%
                     pull())
        } else {
            sort(seq_data_rev() %>%
                     distinct(virus_clean) %>%
                     pull())
        }
        
    })
    
    seq_var_host_genus_rev <- reactive({
        sort(unique(review_data$sequences$host_genus))
    })
    
    seq_var_host_species_rev <- reactive({
        
        if (!is.null(input$host_genus_rev)) {
            sort(seq_data_rev() %>%
                     filter(host_genus %in% input$seq_host_genus_rev) %>%
                     distinct(species) %>%
                     pull())
        } else {
            sort(seq_data_rev() %>%
                     distinct(species) %>%
                     pull())
        }
        
    })
    
    # Filtered data
    seq_data_filtered_rev <- reactive({
        seq_filtered_data <- seq_data_rev()
        
        # Filter based on pathogen family
        if (!is.null(input$seq_pathogen_family_rev) & is.null(input$seq_pathogen_species_rev)) {
            seq_filtered_data <- filter(seq_filtered_data, family %in% input$seq_pathogen_family_rev)
        }
        
        # Filter based on pathogen species
        if (!is.null(input$seq_pathogen_species_rev)) {
            seq_filtered_data <- filter(seq_filtered_data, virus_clean %in% input$seq_pathogen_species_rev)
        }
        
        # Filter based on host genus
        if (!is.null(input$seq_host_genus_rev) & is.null(input$seq_host_species_rev)) {
            seq_filtered_data <- filter(seq_filtered_data, host_genus %in% input$seq_host_genus_rev)
        }
        
        # Filter based on host species
        if (!is.null(input$seq_host_species_rev)) {
            seq_filtered_data <- filter(seq_filtered_data, species %in% input$seq_host_species_rev)
        }
        
        return(seq_filtered_data)
    })
    
    output$seq_review_table_subset <- renderDataTable({
        datatable(seq_data_filtered_rev(), options = list(scrollX = TRUE, autoWidth = TRUE),
                  # colnames = c("Sequence ID", "Rodent ID", "Path ID", "Sampling date", "Study ID", "Host genus", "Host species",
                  #              "Sequence type", "Family", "Virus name", "Coordinate resolution",
                  #              "Lat", "Long", "Accession Number", "Method", "Note", "Date sampled", "Sample location"),
                  class = 'custom-datatable',
                  rownames = FALSE,
                  escape = FALSE)
    })
    
    output$download_seq_data <- downloadHandler(
        filename = function() {
            paste("filtered_seq_data_and_citations_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file) {
            # Create a temporary directory to store files
            temp_dir <- tempdir()
            
            # Filtered table
            write.csv(seq_data_filtered_rev(), file.path(temp_dir, "filtered_seq_data.csv"))
            
            # Filtered citations
            selected_study_ids <- unique(seq_data_filtered_rev()$study_id)
            filtered_citations <- filterCitations(selected_study_ids)
            bib_path <- file.path(temp_dir, "filtered_seq_citations.bib")
            createBibTeXFile(filtered_citations, bib_path)
            
            # Create a zip file with both CSV and BibTeX files
            zip(file, files = c(file.path(temp_dir, "filtered_seq_data.csv"), bib_path), flags = "-j")
        }
    )    

# Plot sequence records ---------------------------------------------------
    
    # Initialize an empty map
    output$sequence_map <- renderLeaflet({
        leaflet() %>%
            addTiles()
    })
    
    colour_palette <- c("darkgreen", "darkred")
    
    # ReactiveVal to track filter status
    filter_applied <- reactiveVal(FALSE)
    
    # Update the map based on filtered data
    observe({
        filtered_data <- seq_data_filtered_rev()
        
        # Check if there are filters applied
        if ((!is.null(input$seq_pathogen_family_rev) && length(input$seq_pathogen_family_rev) > 0) ||
            (!is.null(input$seq_pathogen_species_rev) && length(input$seq_pathogen_species_rev) > 0) ||
            (!is.null(input$seq_host_genus_rev) && length(input$seq_host_genus_rev) > 0) ||
            (!is.null(input$seq_host_species_rev) && length(input$seq_host_species_rev) > 0)) {
            
            # Filter the data based on selected filters
            filtered_data <- filtered_data %>%
                filter(
                    if (!is.null(input$seq_pathogen_family_rev)) family %in% input$seq_pathogen_family_rev else TRUE,
                    if (!is.null(input$seq_pathogen_species_rev)) virus_clean %in% input$seq_pathogen_species_rev else TRUE,
                    if (!is.null(input$seq_host_genus_rev)) host_genus %in% input$seq_host_genus_rev else TRUE,
                    if (!is.null(input$seq_host_species_rev)) associatedTaxa %in% input$seq_host_species_rev else TRUE
                )
            
            # Check if any records are left after filtering
            if (nrow(filtered_data) > 0) {
                # Map organismQuantity to colour
                colours <- ifelse(filtered_data$sequenceType == "Pathogen", "darkred", "darkgreen")
                
                # Clear existing markers and add new ones
                leafletProxy("sequence_map") %>%
                    clearMarkers() %>%
                    addCircleMarkers(data = filtered_data,
                                     lng = ~decimalLongitude,
                                     lat = ~decimalLatitude,
                                     radius = 8,
                                     color = colours,
                                     fillOpacity = 0.7,
                                     popup = ~paste("Study ID: ", study_id, "<br>",
                                                    "Sequence ID: ", sequence_record_id, "<br>",
                                                    "Rodent ID: ", associated_rodent_record_id, "<br>",
                                                    "Pathogen ID: ", associated_pathogen_record_id, "<br>",
                                                    "Accession number: ", accession_number, "<br>",
                                                    "Sequence type: ", sequenceType, "<br>",
                                                    "Rodent: ", species, "<br>",
                                                    "Pathogen: ", virus_clean, "<br>",
                                                    "Coordinate resolution: ", coordinate_resolution, "<br>",
                                                    "Sample date: ", eventDate, "<br>"),
                                     clusterOptions = markerClusterOptions(
                                         spiderfyDistanceMultiplier = 2
                                     )) %>%
                    fitBounds(
                        lat1 = min(filtered_data$decimalLatitude, na.rm = TRUE),
                        lng1 = min(filtered_data$decimalLongitude, na.rm = TRUE),
                        lat2 = max(filtered_data$decimalLatitude, na.rm = TRUE),
                        lng2 = max(filtered_data$decimalLongitude, na.rm = TRUE)
                    )
                
                # Update filter status
                filter_applied(TRUE)
            } else {
                # If no records match the filters, clear markers
                leafletProxy("sequence_map") %>%
                    clearMarkers()
                
                # Update filter status
                filter_applied(FALSE)
            }
        } else {
            # If no filters applied, clear markers
            leafletProxy("sequence_map") %>%
                clearMarkers()
            
            # Update filter status
            filter_applied(FALSE)
        }
    })
    
    
    # Clear markers when filter is removed
    observe({
        if (!filter_applied()) {
            leafletProxy("sequence_map") %>%
                clearMarkers()
        }
    })
    
}


# Run application ---------------------------------------------------------

# Run the application 
shinyApp(ui = ui, server = server)
