# Targeted Data Extraction

# Global parameters and paths

<# 
# Fission Yeast Data
$data_path = "data\GproDIA\fissionyeast"
$library_name = "fissionyeast"
$library_pqp_file = "$library_name.PQP"
$result_name = $library_name
$swath_window_file = "..\swath_window.tsv"
$enable_glycoform_inference = $false
$threads = 4
$mz_extraction_window_ms1 = 10
$mz_extraction_window = 20
#>

# Serum Data
$data_path = "data\GproDIA\serum"
$library_name = "serum"
$library_pqp_file = "$($library_name)_uis.PQP"
$swath_window_file = "..\swath_window.tsv"
$result_name = $library_name
$enable_glycoform_inference = $true
$threads = 4
$mz_extraction_window_ms1 = 10
$mz_extraction_window = 20

# Run OpenSWATH workflow
mkdir "$data_path\result\$result_name" -f | Out-Null

cd "$data_path\mzML"

if (-not $enable_glycoform_inference) {
    ls *.mzML -Name | % {
        OpenSwathWorkflow `
        -Library:retentionTimeInterpretation seconds `
        -RTNormalization:alignmentMethod lowess `
        -mz_extraction_window_ms1 $mz_extraction_window_ms1 -mz_extraction_window_ms1_unit ppm `
        -mz_extraction_window $mz_extraction_window -mz_extraction_window_unit ppm `
        -threads $threads `
        -swath_windows_file "$data_path\$swath_window_file" `
        -tr "..\library\$library_name\$library_pqp_file" `
        -tr_irt "..\library\$library_name\$($library_name)_rtanchor.traML" `
        -in $_ `
        -out_osw "..\result\$result_name\$($_.Replace('.mzML', '.osw'))"
    }
}
else {
    ls *.mzML -Name | % {
        OpenSwathWorkflow `
        -Library:retentionTimeInterpretation seconds `
        -RTNormalization:alignmentMethod lowess `
        -RTNormalization:estimateBestPeptides `
        -mz_extraction_window_ms1 $mz_extraction_window_ms1 -mz_extraction_window_ms1_unit ppm `
        -mz_extraction_window $mz_extraction_window -mz_extraction_window_unit ppm `
        -enable_ms1 true -enable_ipf true `
        -threads $threads `
        -swath_windows_file "$data_path\$swath_window_file" `
        -tr "..\library\$library_name\$library_pqp_file" `
        -tr_irt "..\library\$library_name\$($library_name)_rtanchor.traML" `
        -in $_ `
        -out_osw "..\result\$result_name\$($_.Replace('.mzML', '.osw'))"
    }
}

