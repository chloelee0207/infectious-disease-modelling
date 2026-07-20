# One-off: pull Caldas Novas out of the 129 MB CHIRPS workbook into a small CSV so
# nothing downstream has to touch the full file again. Reads only the columns needed
# (skips the r1h/r3h blocks) and filters to the municipality.
suppressMessages(library(readxl))

# Column order in the sheet:
# date|adm_level|adm_id|PCODE|n_pixels|rfh|rfh_avg|r1h|r1h_avg|r3h|r3h_avg|rfq|r1q|r3q|version|state|municipality
ct <- c("guess","numeric","skip","skip","skip","numeric","numeric",
        "skip","skip","skip","skip","numeric","skip","skip","skip","text","text")

d <- read_excel("rainfall.xlsx", sheet = "bra-rainfall-subnat-5ytd", col_types = ct)
cat("full data sheet rows:", nrow(d), "\n")

ca <- d[!is.na(d$municipality) & tolower(d$municipality) == "caldas novas", ]
cat("Caldas Novas rows:", nrow(ca), "\n")
cat("date range:", format(min(ca$date)), "->", format(max(ca$date)), "\n")
cat("states seen:", paste(unique(ca$state), collapse=", "), "\n")

write.csv(ca, "rainfall_caldas.csv", row.names = FALSE)
cat("Wrote rainfall_caldas.csv\n")
print(utils::head(as.data.frame(ca)))
print(utils::tail(as.data.frame(ca)))
