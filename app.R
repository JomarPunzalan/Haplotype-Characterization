source("run_haplotype.R", local = TRUE)
source("group_review_map.R", local = TRUE)
source("group_review_shiny.R", local = TRUE)

options(shiny.maxRequestSize = max(100 * 1024^2, getOption("shiny.maxRequestSize", 0)))

analysis_defaults <- as.list(formals(run_haplotype))[
  c("max_depth", "split_rule", "min_maf", "fallback_rank",
    "max_missing", "max_heterozygous", "report_max_heterozygous", "marker_pattern")
]

ui <- shiny::fluidPage(
  title = "Haplotype Characterization",
  shiny::tags$head(shiny::tags$style(shiny::HTML("
    body { background: #f5f3f8; color: #231b32; }
    .haplotype-page { max-width: 1180px; margin: 36px auto; padding: 0 12px; }
    .haplotype-page h1 { color: #461d7c; font-weight: 700; margin-bottom: 22px; }
    .haplotype-layout {
      display: grid; grid-template-columns: minmax(275px, 345px) minmax(0, 1fr);
      gap: 20px; align-items: start;
    }
    .haplotype-panel {
      background: white; border: 1px solid #ded5e9; border-top: 4px solid #461d7c;
      border-radius: 12px; padding: 20px; box-shadow: 0 12px 32px #2815420c;
    }
    .haplotype-panel .btn-primary { background: #461d7c; border-color: #461d7c; }
    .haplotype-panel .btn-primary:hover { background: #35145e; border-color: #35145e; }
    .haplotype-panel .btn-preview { color: #185c54; border-color: #27786e; background: #eff9f7; }
    .haplotype-panel .btn-preview:hover { background: #dfefe9; }
    .haplotype-controls h2 { margin: 0 0 12px; }
    .haplotype-controls section { border-top: 1px solid #ded5e9; margin-top: 18px; padding-top: 16px; }
    .haplotype-controls .action-button { width: 100%; white-space: normal; }
    .haplotype-results .tab-pane { padding-top: 16px; }
    .haplotype-downloads { display: flex; flex-wrap: wrap; gap: 8px; margin: 14px 0; }
    .haplotype-summary { margin-top: 22px; }
    .haplotype-summary table { max-width: 420px; }
    .haplotype-qc-table { overflow-x: auto; }
    .review-controls { display: flex; flex-wrap: wrap; gap: 8px 12px; align-items: baseline; }
    .review-controls .form-group { margin-bottom: 8px; }
    .review-inspector { border-top: 1px solid #ded5e9; margin-top: 20px; padding-top: 16px; }
    .review-allele {
      display: inline-block; min-width: 1.8em; border-radius: 4px; padding: 2px 5px;
      text-align: center; font-weight: 600; background: #ece9ef; color: #231b32;
    }
    .review-allele-A { color: #714800; background: #ffeb9c; }
    .review-allele-C { color: #00522a; background: #c6efce; }
    .review-allele-G { color: #002f69; background: #b8cce4; }
    .review-allele-T { color: #861326; background: #ffc7ce; }
    .review-map-legend { display: flex; flex-wrap: wrap; gap: 5px 12px; }
    .review-map-swatch { display: inline-block; width: 0.8em; height: 0.8em; border-radius: 50%; }
    @media (max-width: 770px) {
      .haplotype-page { margin: 22px auto; }
      .haplotype-layout { grid-template-columns: minmax(0, 1fr); }
      .haplotype-panel { padding: 16px; }
    }
  "))),
  shiny::tags$main(
    class = "haplotype-page",
    shiny::tags$h1("Haplotype Characterization"),
    shiny::tags$div(
      class = "haplotype-layout",
      shiny::tags$aside(
        class = "haplotype-panel haplotype-controls",
        shiny::fileInput("snp_file", "SNP file", accept = c(".csv", ".xlsx"), width = "100%"),
        shiny::helpText("Use germplasmName and SNP_chr marker columns. Older files with Line are accepted."),
        shiny::uiOutput("sheet_picker"),
        shiny::tags$section(
          shiny::tags$h2(class = "h4", "SNP filtering"),
          shiny::numericInput(
            "max_missing", "Maximum missing calls (%)", 100 * analysis_defaults$max_missing,
            min = 0, max = 100, step = 1, width = "100%"
          ),
          shiny::numericInput(
            "max_heterozygous", "Maximum unresolved or heterozygous calls for splitting (%)",
            100 * analysis_defaults$max_heterozygous, min = 0, max = 100, step = 1, width = "100%"
          ),
          shiny::numericInput(
            "report_max_heterozygous", "Maximum unresolved or heterozygous calls in filtered report (%)",
            100 * analysis_defaults$report_max_heterozygous, min = 0, max = 100, step = 1, width = "100%"
          ),
          shiny::actionButton("preview_snp_filters", "Preview SNP Filters", class = "btn-preview")
        ),
        shiny::tags$section(
          shiny::tags$h2(class = "h4", "Haplotype analysis"),
          shiny::numericInput(
            "max_depth", "Maximum split levels", analysis_defaults$max_depth, min = 1, step = 1,
            width = "100%"
          ),
          shiny::selectInput(
            "split_rule", "SNP selection rule",
            choices = c(
              "1. Highest MAF" = "highest_maf",
              "2a. Minimum MAF" = "min_maf",
              "2b. Minimum MAF with fallback" = "min_maf_fallback"
            ),
            selected = analysis_defaults$split_rule, width = "100%"
          ),
          shiny::conditionalPanel(
            "input.split_rule === 'min_maf' || input.split_rule === 'min_maf_fallback'",
            shiny::numericInput(
              "min_maf", "Minimum MAF (%)", 100 * analysis_defaults$min_maf, min = 0, max = 50,
              step = 1, width = "100%"
            ),
            shiny::conditionalPanel(
              "input.split_rule === 'min_maf_fallback'",
              shiny::numericInput(
                "fallback_rank", "Top MAF markers to consider", analysis_defaults$fallback_rank,
                min = 1, step = 1, width = "100%"
              )
            ),
            shiny::helpText("Among markers meeting the cutoff, prefer the fewest known homozygous calls.")
          ),
          shiny::actionButton("run_analysis", "Run Haplotype Analysis", class = "btn-primary"),
          shiny::helpText("Uploads run automatically. Use this button to apply changes.")
        )
      ),
      shiny::tags$section(
        class = "haplotype-panel haplotype-results",
        shiny::tabsetPanel(
          id = "results_view",
          shiny::tabPanel(
            "Haplotype Results", value = "haplotype",
            shiny::uiOutput("run_status"),
            shiny::uiOutput("downloads"),
            shiny::uiOutput("group_summary")
          ),
          shiny::tabPanel(
            "SNP Filtering", value = "snp_filtering",
            shiny::uiOutput("qc_status"),
            shiny::uiOutput("qc_download"),
            shiny::tags$div(class = "haplotype-qc-table", shiny::tableOutput("qc_table"))
          ),
          shiny::tabPanel(
            "Group Review", value = "group_review",
            group_review_ui("group_review")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  session_dir <- tempfile("haplotype-session-")
  if (!dir.create(session_dir)) stop("Could not create a temporary directory for this session.")
  session$onSessionEnded(function() unlink(session_dir, recursive = TRUE))

  upload <- shiny::reactiveVal(NULL)
  report <- shiny::reactiveVal(NULL)
  qc_preview <- shiny::reactiveVal(NULL)

  clear_report <- function() {
    previous <- report()
    if (!is.null(previous$directory)) unlink(previous$directory, recursive = TRUE)
    report(NULL)
  }

  shiny::observeEvent(input$snp_file, {
    previous <- upload()
    if (!is.null(previous$directory)) unlink(previous$directory, recursive = TRUE)
    clear_report()
    qc_preview(NULL)

    file <- input$snp_file
    name <- basename(file$name)
    extension <- tolower(tools::file_ext(name))
    if (!extension %in% c("csv", "xlsx")) {
      upload(list(error = "Upload a CSV or XLSX file."))
      return()
    }

    directory <- tempfile("snp-", tmpdir = session_dir)
    path <- file.path(directory, name)
    info <- tryCatch({
      if (!dir.create(directory) || !file.copy(file$datapath, path)) {
        stop("Could not read the uploaded file.", call. = FALSE)
      }
      sheets <- if (extension == "xlsx") haplotype_input_sheets(path) else character()
      if (extension == "xlsx" && !length(sheets)) {
        stop("The workbook needs a worksheet with germplasmName (or legacy Line) and SNP_chr markers.",
             call. = FALSE)
      }
      list(path = path, directory = directory, name = name, sheets = sheets)
    }, error = function(e) {
      unlink(directory, recursive = TRUE)
      list(error = conditionMessage(e))
    })
    upload(info)
  }, priority = 10)

  selected_sheet <- shiny::reactive({
    sheets <- upload()$sheets
    selected <- intersect(input$sheet, sheets)
    if (length(selected)) selected[1L] else if (length(sheets)) sheets[1L]
  })

  output$sheet_picker <- shiny::renderUI({
    sheets <- upload()$sheets
    if (length(sheets) > 1L) {
      shiny::selectInput(
        "sheet", "SNP worksheet", choices = sheets,
        selected = shiny::isolate(selected_sheet())
      )
    }
  })

  number <- function(key, divisor = 1) {
    value <- input[[key]]
    if (is.null(value)) as.numeric(analysis_defaults[[key]]) else as.numeric(value) / divisor
  }

  snp_filters <- shiny::reactive({
    list(
      max_missing = number("max_missing", 100),
      max_heterozygous = number("max_heterozygous", 100),
      report_max_heterozygous = number("report_max_heterozygous", 100)
    )
  })

  analysis_settings <- shiny::reactive({
    rule <- if (is.null(input$split_rule)) analysis_defaults$split_rule else input$split_rule
    c(list(
      max_depth = number("max_depth"),
      split_rule = rule,
      min_maf = if (rule == "highest_maf") analysis_defaults$min_maf else number("min_maf", 100),
      fallback_rank = if (rule == "min_maf_fallback") {
        number("fallback_rank")
      } else {
        analysis_defaults$fallback_rank
      }
    ), snp_filters())
  })

  qc_key <- shiny::reactive({
    list(
      path = upload()$path,
      sheet = selected_sheet(),
      filters = snp_filters()
    )
  })

  shiny::observeEvent(input$preview_snp_filters, {
    shiny::updateTabsetPanel(session, "results_view", selected = "snp_filtering")
    key <- qc_key()
    if (is.null(key$path)) return()
    preview <- tryCatch({
      limits <- unlist(key$filters, use.names = FALSE)
      if (any(!is.finite(limits) | limits < 0 | limits > 1)) {
        stop("Enter SNP filter percentages from 0 to 100.", call. = FALSE)
      }
      snps <- read_haplotype_input(key$path, analysis_defaults$marker_pattern, sheet = key$sheet)
      qc <- do.call(haplotype_marker_qc, c(
        list(calls = snps$calls, unrecognized = snps$unrecognized), key$filters
      ))
      reasons <- cbind(
        "No homozygous calls" = qc$Called == 0,
        "Missing calls above limit" = qc$Missing_Frequency > key$filters$max_missing,
        "Unresolved or heterozygous calls above limit" = qc$Het_Frequency > key$filters$max_heterozygous
      )
      qc$Split_QC <- apply(reasons, 1, function(row) {
        if (any(row)) paste(colnames(reasons)[row], collapse = ", ") else "Passed"
      })
      list(key = key, qc = qc, lines = nrow(snps$data))
    }, error = function(e) {
      list(key = key, error = sub(key$path, upload()$name, conditionMessage(e), fixed = TRUE))
    })
    qc_preview(preview)
  })

  current_qc <- shiny::reactive({
    preview <- qc_preview()
    if (!is.null(preview) && identical(preview$key, qc_key())) preview
  })

  output$qc_status <- shiny::renderUI({
    info <- upload()
    if (!is.null(info$error)) {
      return(shiny::tags$p(role = "alert", class = "text-danger", info$error))
    }
    if (is.null(info$path)) return(shiny::helpText("Upload a SNP file to preview filtering."))
    preview <- current_qc()
    if (is.null(preview)) {
      note <- if (is.null(qc_preview())) {
        "Click Preview SNP Filters to see which markers pass."
      } else {
        "The input or SNP filters changed. Click Preview SNP Filters to refresh."
      }
      return(shiny::helpText(note))
    }
    if (!is.null(preview$error)) {
      return(shiny::tags$p(role = "alert", class = "text-danger", preview$error))
    }
    qc <- preview$qc
    shiny::tags$p(
      role = "status",
      sprintf("%d germplasm entries. %d of %d markers passed splitting QC. %d passed the filtered report criteria.",
              preview$lines, sum(qc$Used_For_Splitting), nrow(qc),
              sum(qc$Included_In_Filtered_Report))
    )
  })

  output$qc_download <- shiny::renderUI({
    if (is.null(current_qc()$qc)) return(NULL)
    shiny::tagList(
      shiny::selectInput(
        "qc_view", "Show markers",
        choices = c("All" = "all", "Excluded from splitting" = "excluded", "Passed splitting QC" = "passed"),
        selected = shiny::isolate(input$qc_view)
      ),
      shiny::downloadButton("download_qc", "Download SNP QC CSV"),
      shiny::helpText("The preview shows up to 40 rows. The CSV includes every marker.")
    )
  })

  output$qc_table <- shiny::renderTable({
    qc <- shiny::req(current_qc()$qc)
    if (identical(input$qc_view, "excluded")) qc <- qc[!qc$Used_For_Splitting, , drop = FALSE]
    if (identical(input$qc_view, "passed")) qc <- qc[qc$Used_For_Splitting, , drop = FALSE]
    shiny::validate(shiny::need(nrow(qc) > 0L, "No markers match this view."))
    qc <- utils::head(qc, 40L)
    data.frame(
      Marker = qc$Marker,
      "Missing (%)" = qc$Missing_Frequency * 100,
      "Unresolved or heterozygous (%)" = qc$Het_Frequency * 100,
      "Split QC" = qc$Split_QC,
      "Filtered report" = ifelse(qc$Included_In_Filtered_Report, "Included", "Excluded"),
      check.names = FALSE
    )
  }, rownames = FALSE)

  output$download_qc <- shiny::downloadHandler(
    filename = function() "snp_filter_preview.csv",
    content = function(file) utils::write.csv(shiny::req(current_qc()$qc), file, row.names = FALSE)
  )

  shiny::observeEvent(list(upload(), selected_sheet(), input$run_analysis), {
    shiny::updateTabsetPanel(session, "results_view", selected = "haplotype")
    info <- upload()
    if (is.null(info$path)) return()

    sheet <- selected_sheet()
    settings <- analysis_settings()
    key <- list(path = info$path, sheet = sheet, settings = settings)
    if (identical(report()$key, key) && is.null(report()$error)) return()
    clear_report()

    directory <- tempfile("report-", tmpdir = session_dir)
    processed <- shiny::withProgress(message = "Characterizing haplotypes", value = 0.2, {
      tryCatch({
        result <- do.call(run_haplotype, c(list(
          input_csv = info$path, output_dir = directory, sheet = sheet, return_review = TRUE
        ), settings))
        paths <- result$paths
        summary <- utils::read.csv(paths[["representatives"]], check.names = FALSE)
        list(
          key = key, directory = directory, paths = paths, sheet = sheet,
          groups = data.frame(Group = summary$Concatenated_Grouping, Lines = summary$Grouping_Count),
          review = result$review, assignments = result$assignments
        )
      }, error = function(e) {
        unlink(directory, recursive = TRUE)
        message <- sub(info$path, info$name, conditionMessage(e), fixed = TRUE)
        list(key = key, error = message)
      })
    })
    report(processed)
  })

  settings_pending <- shiny::reactive({
    result <- report()
    !is.null(result) && !identical(result$key$settings, analysis_settings())
  })

  current_report <- shiny::reactive({
    if (!settings_pending()) report()
  })
  group_review_server("group_review", current_report, settings_pending, report)

  output$run_status <- shiny::renderUI({
    info <- upload()
    if (settings_pending()) {
      return(shiny::tags$p(role = "status", "Settings changed. Click Run Haplotype Analysis to apply them."))
    }
    result <- current_report()
    error <- if (!is.null(info$error)) info$error else result$error
    if (!is.null(error)) {
      return(shiny::tags$p(class = "text-danger", role = "alert", error))
    }
    if (is.null(result$paths)) return(shiny::helpText("Upload a SNP file to generate the haplotype report."))

    source <- if (is.null(result$sheet)) info$name else paste(info$name, result$sheet, sep = " · ")
    lines <- sum(result$groups$Lines)
    groups <- nrow(result$groups)
    shiny::tagList(
      shiny::tags$p(
        role = "status", shiny::tags$strong("Report ready. "),
        sprintf(
          "%s in %s.", ngettext(lines, paste(lines, "line"), paste(lines, "lines")),
          ngettext(groups, paste(groups, "group"), paste(groups, "groups"))
        )
      ),
      shiny::helpText(source),
      if (groups == 1L) shiny::helpText("No groups split under these settings.")
    )
  })

  output$downloads <- shiny::renderUI({
    if (is.null(current_report()$paths)) return(NULL)
    shiny::tags$div(
      class = "haplotype-downloads",
      shiny::downloadButton("download_workbook", "Download Excel report", class = "btn-primary"),
      shiny::downloadButton("download_csvs", "Download CSVs")
    )
  })

  output$group_summary <- shiny::renderUI({
    if (is.null(current_report()$paths)) return(NULL)
    shiny::tags$div(
      class = "haplotype-summary",
      shiny::tags$h2(class = "h4", "Haplotype groups"),
      shiny::tableOutput("groups")
    )
  })

  output$groups <- shiny::renderTable(current_report()$groups, rownames = FALSE)

  output$download_workbook <- shiny::downloadHandler(
    filename = function() "haplotype_report.xlsx",
    content = function(file) {
      path <- shiny::req(current_report()$paths[["workbook"]])
      if (!file.copy(path, file, overwrite = TRUE)) stop("The workbook could not be downloaded.")
    }
  )

  output$download_csvs <- shiny::downloadHandler(
    filename = function() "haplotype_csvs.zip",
    content = function(file) {
      result <- current_report()
      shiny::req(result$paths)
      csvs <- result$paths[setdiff(names(result$paths), "workbook")]
      zip::zipr(file, basename(csvs), root = result$directory)
    }
  )
}

shiny::shinyApp(ui, server)
