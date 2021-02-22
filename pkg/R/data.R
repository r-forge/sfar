#' Spanish Dairy Farm Production
#'
#' @description
#' This dataset contains six years of observations on 247 dairy farms in
#' northern Spain, drawn from 1993-1998. The original data consists of the farm
#' and year identification, plus measurements on one output (i.e. milk), and
#' four inputs (i.e. cows, land, labor and feed).
#'
#' @format A data frame with 1482 observations on the following 29 variables:
#' \describe{
#' \item{FARM}{Farm identification.}
#' \item{AGEL}{Age of the farmer.}
#' \item{YEAR}{Year identification.}
#' \item{COWS}{Number of milking cows.}
#' \item{LAND}{Agricultural area.}
#' \item{MILK}{Milk production.}
#' \item{LABOR}{Labor.}
#' \item{FEED}{Feed.}
#' \item{YIT}{Log of \code{MILK}.}
#' \item{X1}{Log of \code{COWS}.}
#' \item{X2}{Log of \code{LAND}.}
#' \item{X3}{Log of \code{LABOR}.}
#' \item{X4}{Log of \code{FEED}.}
#' \item{X11}{1/2 * \code{X1}^2.}
#' \item{X22}{1/2 * \code{X2}^2.}
#' \item{X33}{1/2 * \code{X3}^2.}
#' \item{X44}{1/2 * \code{X4}^2.}
#' \item{X12}{\code{X1} * \code{X2}.}
#' \item{X13}{\code{X1} * \code{X3}.}
#' \item{X14}{\code{X1} * \code{X4}.}
#' \item{X23}{\code{X2} * \code{X3}.}
#' \item{X24}{\code{X2} * \code{X4}.}
#' \item{X34}{\code{X3} * \code{X4}.}
#' \item{YEAR93}{Dummy for \code{YEAR = 1993}.}
#' \item{YEAR94}{Dummy for \code{YEAR = 1994}.}
#' \item{YEAR95}{Dummy for \code{YEAR = 1995}.}
#' \item{YEAR96}{Dummy for \code{YEAR = 1996}.}
#' \item{YEAR97}{Dummy for \code{YEAR = 1997}.}
#' \item{YEAR98}{Dummy for \code{YEAR = 1998}.}
#' }
#' @details This dataset has been used in Alvarez \emph{et al.} (2004). The data
#'   have been normalized so that the logs of the inputs sum to zero over the
#'   1482 observations.
#' @source
#'   \url{http://pages.stern.nyu.edu/~wgreene/Econometrics/oldPanelDataSets.htm}
#' @references Alvarez, A., C. Arias, and W. Greene. 2004. Accounting for
#'   unobservables in production models: management and inefficiency.
#'   \emph{Econometric Society}, 1--20.
#' @examples str(dairyspain)
#' summary(dairyspain)
#' @keywords datasets
"dairyspain"


#' U.S. electric power generation
#'
#' This dataset is on U.S. electic power generation.
#'
#' @format A data frame with 123 observations on the following 9 variables:
#'   \describe{
#'   \item{firm}{Firm identification.}
#'   \item{cost}{Total cost in 1970, MM USD.}
#'   \item{output}{Output in million KwH.}
#'   \item{lprice}{Labor price.}
#'   \item{lshare}{Labor's cost share.}
#'   \item{cprice}{Capital price.}
#'   \item{cshare}{Capital's cost share.}
#'   \item{fprice}{Fuel price.}
#'   \item{fshare}{Fuel's cost share.}
#'   }
#' @details The data set is from Christensen and Greene (1976) and has also been
#'   used in Greene (1990).
#' @source \url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
#' @references
#' Christensen, L.R., and W.H. Greene. 1976. Economies of scale in US electric
#' power generation. \emph{The Journal of Political Economy}, 655--676.
#'
#' Greene, W.H. 1990. A Gamma-Distributed Stochastic Frontier Model.
#' \emph{Journal of econometrics}, \bold{46}:141--163.
#' @examples str(electricity)
#' summary(electricity)
#' @keywords datasets
"electricity"


#' Rice Production in the Philippines
#'
#' This dataset contains annual data collected from 43 smallholder rice
#' producers in the Tarlac region of the Philippines between 1990 and 1997.
#'
#' @format A data frame with 344 observations on the following 17 variables:
#' \describe{
#' \item{YEARDUM}{Time period (1= 1990, ..., 8 = 1997).}
#' \item{FMERCODE}{Farmer code (1, ..., 43).}
#' \item{PROD}{Output (tonnes of freshly threshed rice).}
#' \item{AREA}{Area planted (hectares).}
#' \item{LABOR}{Labour used (man-days of family and hired labour).}
#' \item{NPK}{Fertiliser used (kg of active ingredients).}
#' \item{OTHER}{Other inputs used (Laspeyres index = 100 for Firm 17 in 1991).}
#' \item{PRICE}{Output price (pesos per kg).}
#' \item{AREAP}{Rental price of land (pesos per hectare).}
#' \item{LABORP}{Labour price (pesos per hired man-day.}
#' \item{NPKP}{Fertiliser price (pesos per kg of active ingredient).}
#' \item{OTHERP}{Price of other inputs (implicit price index).}
#' \item{AGE}{Age of the household head (years).}
#' \item{EDYRS}{Education of the household head (years).}
#' \item{HHSIZE}{Household size.}
#' \item{NADULT}{Number of adults in the household.}
#' \item{BANRAT}{Percentage of area classified as bantog (upland) fields.}
#' }
#' @details This data set is published as supplement to Coelli \emph{et al.}
#'   (2005). While most variables of this data set were supplied by the
#'   International Rice Research Institute (IRRI), some were calculated by
#'   Coelli \emph{et al.} (2005, see p. 325--326). The survey is described in
#'   Pandey \emph{et al.} (1999).
#' @source Supplementary files for Coelli \emph{et al.} (2005),
#'   \url{http://www.uq.edu.au/economics/cepa/crob2005/software/CROB2005.zip}.
#'   See also \code{\link[frontier]{riceProdPhil}}.
#' @references Coelli, T. J., Rao, D. S. P., O'Donnell, C. J., and Battese, G.
#'   E. (2005). \emph{An Introduction to Efficiency and Productivity Analysis},
#'   Springer, New York.
#'
#'   Pandey, S., Masciat, P., Velasco, L, and Villano, R. (1999). Risk analysis
#'   of a rainfed rice production system system in Tarlac, Central Luzon,
#'   Philippines. \emph{Experimental Agriculture}, \bold{35}:225--237.
#' @examples str(ricephil)
#' summary(ricephil)
#' @keywords datasets
"ricephil"


#' Swiss Railways
#'
#' This dataset is an unbalanced panel of 50 Swiss railway companies over the
#' period 1985-1997.
#'
#' @format A data frame with 605 observations on the following 30 variables:
#' \describe{
#' \item{ID}{Firm identification.}
#' \item{YEAR}{Year identification.}
#' \item{NI}{Number of years observed.}
#' \item{STOPS}{Number of stops in network.}
#' \item{NETWORK}{Network lenght (in meters).}
#' \item{NARROW_T}{Dummy variable for railroads with narrow track.}
#' \item{RACK}{Dummy variable for ‘rack rail’ in network.}
#' \item{TUNNEL}{Dummy variable for network with tunnels over 300 meters on
#' average.}
#' \item{T}{Time indicator, first year = 0.}
#' \item{Q2}{Passenger output – passenger km.}
#' \item{Q3}{Freight output – ton km.}
#' \item{CT}{Total cost (1,000 Swiss franc).}
#' \item{PL}{Labor price.}
#' \item{PE}{Electricity price.}
#' \item{PK}{Capital price.}
#' \item{VIRAGE}{1 for railroads with curvy tracks.}
#' \item{LNCT}{Log of \code{CT}/\code{PE}.}
#' \item{LNQ2}{Log of \code{Q2}.}
#' \item{LNQ3}{Log of \code{Q3}.}
#' \item{LNNET}{Log of \code{NETWORK}/1000.}
#' \item{LNPL}{Log of \code{PL}/\code{PE}.}
#' \item{LNPE}{Log of \code{PE}.}
#' \item{LNPK}{Log of \code{PK}/\code{PE}.}
#' \item{LNSTOP}{Log of \code{STOPS}.}
#' \item{MLNQ2}{Cross sectional mean of \code{LNQ2}.}
#' \item{MLNQ3}{Cross sectional mean of \code{LNQ3}.}
#' \item{MLNNET}{Cross sectional mean of \code{LNNET}.}
#' \item{MLNPL}{Cross sectional mean of \code{LNPL}.}
#' \item{MLNPK}{Cross sectional mean of \code{LNPK}.}
#' \item{MLNSTOP}{Cross sectional mean of \code{LNSTOP}.}
#' }
#' @details The dataset is extracted from the annual reports of the Swiss
#'   Federal Office of Statistics on public transport companies and has been
#'   used in Farsi`\emph{et al.} (2005).
#' @source
#'   \url{http://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm}
#'   \url{http://people.stern.nyu.edu/wgreene/Microeconometrics.htm}
#' @references Farsi, M., M. Filippini, and W. Greene. 2005. Efficiency
#'   Measurement in Network Industries: Application to the Swiss Railway
#'   Companies. \emph{Journal of Regulatory Economics}, \bold{28}:69--90.
#' @examples str(swissrailways)
#' summary(swissrailways)
#' @keywords datasets
"swissrailways"


#' U.S. electricity generating plants
#'
#' This dataset contains data on fossil fuel fired steam electric
#' power generation plants in the United States between 1986 and 1996.
#'
#' @format  A data frame with 791 observations on the following 11 variables:
#' \describe{
#' \item{\code{firm}}{Plant identification.}
#' \item{\code{year}}{Year identification.}
#' \item{\code{y}}{Net-steam electric power generation in megawatt-hours.}
#' \item{\code{regu}}{Dummy variable which takes a value equal to 1 if the power
#' plant is in a state which enacted legislation or issued a regulatory order to
#' implement retail access during the sample period, and a value equal to 0
#' otherwise.}
#' \item{\code{k}}{Capital stock.}
#' \item{\code{labor}}{Labor and maintenance.}
#' \item{\code{fuel}}{Fuel.}
#' \item{\code{wl}}{Labor price.}
#' \item{\code{wf}}{Fuel price.}
#' \item{\code{wk}}{Capital price.}
#' \item{\code{tc}}{Total cost.}
#' }
#' @details The data set has been used in Kumbhakar \emph{et al.} (2014).
#' @source
#'   \url{https://sites.google.com/site/sfbook2014/home/for-stata-v12-v13-v14}
#' @references Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. A
#'   Practitioner's Guide to Stochastic Frontier Analysis Using Stata: Cambridge
#'   University Press.
#' @examples str(utility)
#' summary(utility)
#' @keywords datasets
"utility"


#' World production frontier
#'
#' This dataset provide information on production related variables
#' for eighty-two countries over the period 1960–1987.
#'
#' @format A data frame with 2296 observations on the following 12 variables:
#' \describe{
#' \item{\code{country}}{Country name.}
#' \item{\code{code}}{Country identification.}
#' \item{\code{yr}}{Year identification.}
#' \item{\code{y}}{GDP in 1987 U.S. dollars.}
#' \item{\code{k}}{Physical capital stock in 1987 U.S. dollars.}
#' \item{\code{l}}{Labor (number of individuals in the workforce between the
#' ages of 15 and 64).}
#' \item{\code{h}}{Human capital-adjusted labor.}
#' \item{\code{ly}}{Log of \code{y}.}
#' \item{\code{lk}}{Log of \code{k}.}
#' \item{\code{ll}}{Log of \code{l}.}
#' \item{\code{lh}}{Log of \code{h}.}
#' \item{\code{initStat}}{Log of the initial capital to labor ratio of each
#' country, \code{lk} - \code{ll}, measured at the beginning of the sample
#' period.} }
#' @details The dataset is from the World Bank STARS database and has been used
#'   in Kumbhakar \emph{et al.} (2014).
#' @source
#'   \url{https://sites.google.com/site/sfbook2014/home/for-stata-v12-v13-v14}
#' @references Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. A
#'   Practitioner's Guide to Stochastic Frontier Analysis Using Stata: Cambridge
#'   University Press.
#' @examples str(worldprod)
#' summary(worldprod)
#' @keywords datasets
"worldprod"
