# ---- missing interactions ----
# install.packages("rglobi")
library(rglobi)

# e.g. get all “eats” interactions where the source is "Panthera leo"
ints_euporbia <- get_interactions_by_taxa(
  sourcetaxon   = "euphorbia balsamifera"
)

# Take a look
view(ints_euporbia)

ints_camponotus <- get_interactions_by_taxa(
  sourcetaxon   = "camponotus"
)

# Take a look
view(ints_camponotus)

ints_feae <- get_interactions_by_taxa(
  sourcetaxon   = "feae"
)

# Take a look
view(ints_feae)

ints_camponotus <- ints_feae %>% filter(source_taxon_name == "Camponotus feae")
view(ints_camponotus)

ints_anaspis <- get_interactions_by_taxa(
  sourcetaxon   = "anaspis proteus"
)

# Take a look
view(ints_anaspis)

ints_urell <- get_interactions_by_taxa(
  sourcetaxon   = "urelliosoma"
)

# Take a look
view(ints_urell)

# for Brazil
ints_bombus <- get_interactions_by_taxa(
  sourcetaxon   = "cuphea ericoides"
)

# Take a look
view(ints_bombus)

#after corrections
ints_allaudi <- get_interactions_by_taxa(
  sourcetaxon   = "anthophora alluaudi"
)

# Take a look
view(ints_allaudi)

ints_camponotus <- get_interactions_by_taxa(
  sourcetaxon   = "camponotus feae"
)

# Take a look
view(ints_camponotus)

ints_euphorbia <- get_interactions_by_taxa(
  sourcetaxon   = "euphorbia balsamifera"
)

# Take a look
view(ints_euphorbia)
