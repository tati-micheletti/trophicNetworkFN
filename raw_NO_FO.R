###############################################################################
## 7. Read summarized prey data -- "Newer Data"
###############################################################################

raw_psiri <- fread("data/raw_data_NO_FO.csv")

setnames(
  raw_psiri,
  old = c(
    "Species Eaten",
    "Species Evaluated",
    "FO",
    "NO",
    "Volume"
  ),
  new = c(
    "prey",
    "predator",
    "fo",
    "count",
    "volume"
  )
)

## Remove predator species from the prey list
raw_psiri <- raw_psiri[
  !prey %in% c(
    "Teiu",
    "Gato",
    "Rato",
    "Sapo"
  )
]

## Replace missing values with zero
raw_psiri[is.na(fo),     fo := 0]
raw_psiri[is.na(count),  count := 0]
raw_psiri[is.na(volume), volume := 0]
