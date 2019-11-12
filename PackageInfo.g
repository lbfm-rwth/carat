#############################################################################
##  
##  Demo PackageInfo.g for the GitHubPagesForGAP
##

SetPackageInfo( rec(

PackageName := "Carat",

Subtitle := "Crystallographic AlgoRithms And Tables",
Version := "1.1",
Date := "10/11/2019", # dd/mm/yyyy format
License := "GPL v@.0",

Persons := [
  #TODO: Add original authors
  rec(
    LastName      := "Horn",
    FirstNames    := "Max",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "max.horn@uni-siegen.de",
    WWWHome       := "https://www.quendi.de/math",
    PostalAddress := Concatenation(
                       "Department Mathematik\n",
                       "Universität Siegen\n",
                       "Walter-Flex-Straße 3\n",
                       "57072 Siegen\n",
                       "Germany" ),
    Place         := "Siegen",
    Institution   := "Universität Siegen"
  ),
  rec(
    LastName      := "Berhardt",
    FirstNames    := "Dominik",
    IsAuthor      := false,
    IsMaintainer  := true,
    Email         := "carat@momo.math.rwth-aachen.de",
    WWWHome       := "https://www.mathb.rwth-aachen.de/cms/MATHB/Der-Lehrstuhl/Team/Wissenschaftliche-Beschaeftigte/~rnsg/Dominik-Bernhardt/lidx/1/",
    PostalAddress := Concatenation(
                       "RWTH Aachen University\n",
                       "Lehrstuhl B für Mathematik\n",
                       "Pontdriesch 10-16\n",
                       "52062 Aachen\n",
                       "Germany" ),
    Place         := "Aachen",
    Institution   := "RWTH Aachen University"
  ),
],

Status := "other",

# The following are not strictly necessary in your own PackageInfo.g
# (in the sense that update.g only looks at the usual fields
# like PackageWWWHome, ArchiveURL etc.). But they are convenient
# if you use exactly the scheme for your package website that we propose.
GithubUser := "lbfm-rwth",
GithubRepository := ~.PackageName,
GithubWWW := Concatenation("https://github.com/", ~.GithubUser, "/", ~.GithubRepository),

PackageWWWHome := Concatenation("https://", ~.GithubUser, ".github.io/", ~.GithubRepository, "/"),
README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
# The following assumes you are using the Github releases system. If not, adjust
# it accordingly.
ArchiveURL     := Concatenation(~.GithubWWW,
                    "/releases/download/v", ~.Version, "/",
                    ~.GithubRepository, "-", ~.Version),

ArchiveFormats := ".tar.gz .tar.bz2",

AbstractHTML := 
"Carat is a computer package which\
handles enumeration, construction, recognition,\
and comparison problems for crystallographic groups\
up to dimension 6. The name CARAT itself is an\
acronym for Crystallographic AlgoRithms And Tables.",

#TODO: Put .pdf into gh-pages
PackageDoc := rec(
  BookName  := "GitHubPagesForGAP",
  ArchiveURLSubset := ["tex"],
  HTMLStart := "tex/index.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "A set of C programms to deal with crystallographic groups",
),

#TODO: What are our dependencies? 
# This is just here to make the site work.
Dependencies := rec(
  GAP := ">=4.8.1",
  NeededOtherPackages := [
    ["GAPDoc", ">= 1.2"],
    ["IO", ">= 4.1"],
  ],
  SuggestedOtherPackages := [["orb", ">= 4.2"]],
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,

Keywords := ["GitHub Pages", "CARAT"]
));


