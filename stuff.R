data("frogs_datatable")

frogs_datalist <- DAISIE_dataprep(
  datatable = frogs_datatable,
  island_age = 30,
  M = 300
)

DAISIE::DAISIE_ML(
  datalist = frogs_datalist,
  initparsopt = c(0.18,0.03,0.0006,2),
  idparsopt = c(1,2,4,5),
  ddmodel = 0,
  parsfix = Inf,
  idparsfix = 3
)
