# Using a Google Sheet

The app accepts uploaded CSV or XLSX files, not a Google Sheets URL or
published link. It reads uploaded files only; it does not fetch sheet URLs.

In Google Sheets, choose **File > Download** and export as CSV or Microsoft
Excel (.xlsx). Then upload that downloaded file in the app. A CSV contains one
sheet, so download each sheet you need as a separate CSV; an XLSX workbook can
carry multiple sheets for selection in the app.

Make sure the downloaded table has a header row and map the correct header row
in the app. Map each distinct sex value to male, female, or unknown. For
X-linked genes, enter male genotypes as `A/Y` and female genotypes as `A/a`.
