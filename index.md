---
layout: default
---

# CARAT

Carat is a computer package which handles enumeration, construction,
recognition and comparison problems for crystallographic groups up to dimension 6.
The name CARAT itself is anacronym for Crystallographic AlgoRithms And Tables.

The current version of this package is version {{site.data.release.version}}, released on {{site.data.release.date}}.
For more information, please refer to [**the package manual**](doc/index.html).
There is also a [README](README.html) file.

## Dependencies

The following software is needed:
- A C compiler, for example [gcc](https://gcc.gnu.org)
- The [GMP library](https://gmplib.org)


## Authors

Jürgen Opgenorth,
[Wilhelm Plesken](https://www.mathb.rwth-aachen.de/cms/MATHB/Der-Lehrstuhl/Team/Professorinnen-und-Professoren/~rnmq/Univ-Prof-i-R-Dr-rer-nat-Wilhelm-Pl/lidx/1/),
Tilman Schulz.

## Contributors

[Dominik Bernhardt](https://www.mathb.rwth-aachen.de/cms/MATHB/Der-Lehrstuhl/Team/Wissenschaftliche-Beschaeftigte/~rnsg/Dominik-Bernhardt/lidx/1/),
[Franz Gähler](https://www.math.uni-bielefeld.de/~gaehler/),
[Max Horn](https://www.quendi.de/math).


## Citing

Please, cite this package as

<p class='BibEntry'>
[<span class='BibKey'>OPS19</span>]   <b class='BibAuthor'>Opgenorth, J., Plesken, W. and Schulz, T.</b>,
 <i class='BibTitle'>Carat, Crystallographic AlgoRithms And Tables,
         Version {{site.data.release.version}}</i>
 (<span class='BibYear'>{{site.data.release.date | date: "%Y"}}</span>)<br />
<span class='BibHowpublished'><a href="https://github.com/lbfm-rwth/carat/">https://github.com/lbfm-rwth/carat/</a></span>.
</p>

Here is a BibTeX entry:
<pre id='bibtex'>
@misc{ Carat{{site.data.release.version}},
  author =           {Opgenorth, J., Plesken, W. and Schulz, T.},
  title =            {% raw %}{{% endraw %}{Carat}, Crystallographic AlgoRithms And Tables,
                      {V}ersion {{site.data.release.version}}},
  month =            {% raw %}{{% endraw %}{{site.data.release.date | date: "%b"}}},
  year =             {% raw %}{{% endraw %}{{site.data.release.date | date: "%Y"}}},
  howpublished =     {\url{https://github.com/lbfm-rwth/carat/}},
  printedkey =       {OPS{{site.data.release.date | date: "%y"}}}
}
</pre>

<button type="button" class="btn active" data-clipboard-target="#bibtex">Copy BibTeX to clipboard</button>


## Feedback

For bug reports, feature requests and suggestions, please use the
[issue tracker](https://github.com/lbfm-rwth/carat/issues).
