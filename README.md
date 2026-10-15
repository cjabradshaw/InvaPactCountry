# InvaPact synthesis country-scale analysis
<img align="right" src="https://github.com/InvaPact-Synthesis/.github/blob/main/profile/www/InvaPact%20logo-transp.png" alt="InvaPact logo-transp.png" width="250" style="margin-top: 20px">

Drivers of among-country InvaPact scores

#### Investigators
- <a href="">Laís Carneiro</a>, Université Paris-Saclay
- <a href="https://scholar.google.com/citations?user=dH_qGsMAAAAJ&hl=es">Elena Angulo</a>, Consejo Superior de Investigaciones Científicas
- <a href="https://globalecologyflinders.com/people/#CJAB">Corey Bradshaw</a>, Flinders University
- <a href="https://www.biodiversitydynamics.fr/franck-courchamp/">Franck Courchamp</a>, Université Paris-Saclay
- <a href="https://pure.qub.ac.uk/en/persons/ross-cuthbert-2">Ross Cuthbert</a>, Queen's University Belfast
- <a href="https://scholar.google.com/citations?user=Nyeg3eIAAAAJ&hl=fr">Christophe Diagne</a>, Institut de recherche pour le développement
- <a href="">Daniel Pincheira-Donoso</a>, Queen's University Belfast
- <a href="">Boris Leroy</a>, Muséum d'Histoire naturelle

<br>
Accompanies paper:<br>
<br>
Carneiro L, <em>et al</em>. 2026. Global patterns and drivers of the ecological impacts of invasive species. in preparation

## Scripts
- <code>InvaPactCountry.R</code>: main R code for analysis

## Data
- <em>gbif_records_country_year.csv</em>: <a href="https://www.gbif.org/">Global Biodiversity Information Facility</a> (GBIF) # records per country
- <em>WB_WDI_SI_POV_GINI.csv</em>: latest Gini coefficient (source: <a href="https://databank.worldbank.org">World Bank</a>
- <em>HDIPP.csv</em>: planetary pressure-adjusted Human Development Index (source: <a href="https://hdr.undp.org/data-center/human-development-index#/indicies/HDI">United Nations Development Programme</a>
- <em>IGS.csv</em>: imports of goods and services (source: <a href="https://databank.worldbank.org">World Bank</a>
- <em>cpi24.csv</em>: Corruption Perceptions Index (2024) (source: <a href="https://www.transparency.org/">Transparency International</a>
- <em>NEFBA_Data_global_national_ecological_footprint_biocapacity.xlsx</em>: ecological footprint (source: <a href="http://footprintnetwork.org">Footprint Network</a>
- <em>freedom2025.csv</em>: aggregate freedom index (source: <a href="https://freedomhouse.org/report/freedom-world#Data">Freedom House</a>
- <em>popXage.csv</em>: population size (source: <a href="https://population.un.org/wpp/">United Nations Population Division World Population Prospects</a>
- <em>DCWI.csv</em>: domestic comprehensive wealth index (source: <a href="https://data360.worldbank.org/">World Bank</a>
- <em>gdppcppp.csv</em>: purchasing power parity-adjusted per-capita gross national product (source: <a href="https://data360.worldbank.org/">World Bank</a>
- <em>cont_cntry.csv</em>: country and region codes

## Required R packages
<code>boot</code>, <code>countrycode</code>, <code>dismo</code>, <code>dplyr</code>, <code>gbm</code>, <code>ggplot2</code>, <code>invapact</code>, <code>readxl</code>
