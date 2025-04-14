<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en"><head>

<meta charset="utf-8">
<meta name="generator" content="quarto-1.3.450">

<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">


<title>readme</title>
<style>
code{white-space: pre-wrap;}
span.smallcaps{font-variant: small-caps;}
div.columns{display: flex; gap: min(4vw, 1.5em);}
div.column{flex: auto; overflow-x: auto;}
div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
ul.task-list{list-style: none;}
ul.task-list li input[type="checkbox"] {
  width: 0.8em;
  margin: 0 0.8em 0.2em -1em; /* quarto-specific, see https://github.com/quarto-dev/quarto-cli/issues/4556 */ 
  vertical-align: middle;
}
/* CSS for syntax highlighting */
pre > code.sourceCode { white-space: pre; position: relative; }
pre > code.sourceCode > span { display: inline-block; line-height: 1.25; }
pre > code.sourceCode > span:empty { height: 1.2em; }
.sourceCode { overflow: visible; }
code.sourceCode > span { color: inherit; text-decoration: inherit; }
div.sourceCode { margin: 1em 0; }
pre.sourceCode { margin: 0; }
@media screen {
div.sourceCode { overflow: auto; }
}
@media print {
pre > code.sourceCode { white-space: pre-wrap; }
pre > code.sourceCode > span { text-indent: -5em; padding-left: 5em; }
}
pre.numberSource code
  { counter-reset: source-line 0; }
pre.numberSource code > span
  { position: relative; left: -4em; counter-increment: source-line; }
pre.numberSource code > span > a:first-child::before
  { content: counter(source-line);
    position: relative; left: -1em; text-align: right; vertical-align: baseline;
    border: none; display: inline-block;
    -webkit-touch-callout: none; -webkit-user-select: none;
    -khtml-user-select: none; -moz-user-select: none;
    -ms-user-select: none; user-select: none;
    padding: 0 4px; width: 4em;
  }
pre.numberSource { margin-left: 3em;  padding-left: 4px; }
div.sourceCode
  {   }
@media screen {
pre > code.sourceCode > span > a:first-child::before { text-decoration: underline; }
}
</style>


<script src="README_files/libs/clipboard/clipboard.min.js"></script>
<script src="README_files/libs/quarto-html/quarto.js"></script>
<script src="README_files/libs/quarto-html/popper.min.js"></script>
<script src="README_files/libs/quarto-html/tippy.umd.min.js"></script>
<script src="README_files/libs/quarto-html/anchor.min.js"></script>
<link href="README_files/libs/quarto-html/tippy.css" rel="stylesheet">
<link href="README_files/libs/quarto-html/quarto-syntax-highlighting.css" rel="stylesheet" id="quarto-text-highlighting-styles">
<script src="README_files/libs/bootstrap/bootstrap.min.js"></script>
<link href="README_files/libs/bootstrap/bootstrap-icons.css" rel="stylesheet">
<link href="README_files/libs/bootstrap/bootstrap.min.css" rel="stylesheet" id="quarto-bootstrap" data-mode="light">


</head>

<body class="fullcontent">

<div id="quarto-content" class="page-columns page-rows-contents page-layout-article">

<main class="content" id="quarto-document-content">



<!-- README.md is generated from README.qmd via: quarto render README.qmd --output README.md -->
<section id="gamlss2-infrastructure-for-flexible-distributional-regression" class="level1">
<h1>gamlss2: Infrastructure for Flexible Distributional Regression</h1>
<section id="overview" class="level2">
<h2 class="anchored" data-anchor-id="overview">Overview</h2>
<p>The primary purpose of this package is to facilitate the creation of advanced infrastructures designed to enhance the GAMLSS modeling framework. Notably, the <code>gamlss2</code> package represents a significant overhaul of its predecessor, <a href="https://cran.r-project.org/package=gamlss"><code>gamlss</code></a>, with a key emphasis on improving estimation speed and incorporating more flexible infrastructures. These enhancements enable the seamless integration of various algorithms into GAMLSS, including gradient boosting, Bayesian estimation, regression trees, and forests, fostering a more versatile and powerful modeling environment.</p>
<p>Moreover, the package expands its compatibility by supporting all model terms from the base R <a href="https://cran.r-project.org/package=mgcv"><code>mgcv</code></a> package. Additionally, the <code>gamlss2</code> package introduces the capability to accommodate more than four parameter families. Essentially, this means that users can now specify any type of model using these new infrastructures, making the package highly flexible and accommodating to a wide range of modeling requirements.</p>
<ul>
<li>The main model function is <a href="https://gamlss-dev.github.io/gamlss2/man/gamlss2.html"><code>gamlss2()</code></a>.</li>
<li>The default optimizer functions is <a href="https://gamlss-dev.github.io/gamlss2/man/RS_CG.html"><code>RS()</code></a>. Optimizer functions can be exchanged.</li>
<li>Most important methods: <code>summary()</code>, <a href="https://gamlss-dev.github.io/gamlss2/man/plots.html"><code>plot()</code></a>, <a href="https://gamlss-dev.github.io/gamlss2/man/predict.gamlss2.html"><code>predict()</code></a>.</li>
<li>Easy development of new family objects, see <a href="https://gamlss-dev.github.io/gamlss2/man/gamlss2.family.html"><code>?gamlss2,family</code></a>.</li>
<li>User-specific “special” terms are possible, see <a href="https://gamlss-dev.github.io/gamlss2/man/special_terms.html"><code>?special_terms</code></a>.</li>
</ul>
<p>For examples, please visit the manual pages.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb1"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb1-1"><a href="#cb1-1" aria-hidden="true" tabindex="-1"></a><span class="fu">help</span>(<span class="at">package =</span> <span class="st">"gamlss2"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
</section>
<section id="installation" class="level2">
<h2 class="anchored" data-anchor-id="installation">Installation</h2>
<p>The development version of <code>gamlss2</code> can be installed via</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb2"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a><span class="fu">install.packages</span>(<span class="st">"gamlss2"</span>,</span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a>  <span class="at">repos =</span> <span class="fu">c</span>(<span class="st">"https://gamlss-dev.R-universe.dev"</span>,</span>
<span id="cb2-3"><a href="#cb2-3" aria-hidden="true" tabindex="-1"></a>            <span class="st">"https://cloud.R-project.org"</span>))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
</section>
<section id="licence" class="level2">
<h2 class="anchored" data-anchor-id="licence">Licence</h2>
<p>The package is available under the <a href="https://www.gnu.org/licenses/gpl-3.0.html">General Public License version 3</a> or <a href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.html">version 2</a></p>
</section>
<section id="illustration" class="level2">
<h2 class="anchored" data-anchor-id="illustration">Illustration</h2>
<p>The package is designed to follow the workflow of well-established model fitting functions like <code>lm()</code> or <code>glm()</code>, i.e., the step of estimating full distributional regression models is actually not very difficult.</p>
<p>To illustrate the workflow using <code>gamlss2</code>, we analyze the <code>WeatherGermany</code> data,</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb3"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a><span class="cf">if</span>(<span class="sc">!</span>(<span class="st">"WeatherGermany"</span> <span class="sc">%in%</span> <span class="fu">installed.packages</span>())) {</span>
<span id="cb3-2"><a href="#cb3-2" aria-hidden="true" tabindex="-1"></a>  <span class="fu">install.packages</span>(<span class="st">'WeatherGermany'</span>,</span>
<span id="cb3-3"><a href="#cb3-3" aria-hidden="true" tabindex="-1"></a>    <span class="at">repos =</span> <span class="fu">c</span>(<span class="st">"https://gamlss-dev.r-universe.dev"</span>,</span>
<span id="cb3-4"><a href="#cb3-4" aria-hidden="true" tabindex="-1"></a>              <span class="st">"https://cloud.r-project.org"</span>))</span>
<span id="cb3-5"><a href="#cb3-5" aria-hidden="true" tabindex="-1"></a>}</span>
<span id="cb3-6"><a href="#cb3-6" aria-hidden="true" tabindex="-1"></a><span class="fu">data</span>(<span class="st">"WeatherGermany"</span>, <span class="at">package =</span> <span class="st">"WeatherGermany"</span>)</span>
<span id="cb3-7"><a href="#cb3-7" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(WeatherGermany)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>  id       date Wmax pre Tmax Tmin sun name alt     lat    lon
1  1 1981-01-01   NA 1.7  3.4 -5.0  NA Aach 478 47.8413 8.8493
2  1 1981-01-02   NA 1.7  1.2 -0.4  NA Aach 478 47.8413 8.8493
3  1 1981-01-03   NA 5.4  5.4  1.0  NA Aach 478 47.8413 8.8493
4  1 1981-01-04   NA 8.8  5.6 -0.4  NA Aach 478 47.8413 8.8493
5  1 1981-01-05   NA 3.7  1.2 -2.4  NA Aach 478 47.8413 8.8493
6  1 1981-01-06   NA 4.0  1.2 -2.2  NA Aach 478 47.8413 8.8493</code></pre>
</div>
</div>
<p>The dataset contains daily observations from weather stations across Germany. It includes the station identifier (<code>id</code>), the recording <code>date</code>, the maximum wind speed (<code>Wmax</code>, in m/s), the precipitation amount (<code>pre</code>, in mm), the maximum and minimum temperatures (<code>Tmax</code> and <code>Tmin</code>, in °C), and the number of sunshine hours (<code>sun</code>). Additionally, it provides the station’s <code>name</code>, <code>alt</code>itude (in meters above sea level), and its geographic coordinates (<code>lon</code>gitude and <code>lat</code>itude).</p>
<p>In this example, we use daily maximum temperature (<code>Tmax</code>) data to estimate a climatology model based on over 37 years of observations from Germany’s highest meteorological station, located at Zugspitze. Situated at an altitude of 2956 meters above sea level, this station provides a unique dataset for high-altitude climate analysis.</p>
<p>First, we subset the dataset to include only observations from the Zugspitze station.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb5"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb5-1"><a href="#cb5-1" aria-hidden="true" tabindex="-1"></a>d <span class="ot">&lt;-</span> <span class="fu">subset</span>(WeatherGermany, name <span class="sc">==</span> <span class="st">"Zugspitze"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>Before estimating a climatology model using <code>gamlss2</code>, it is good practice to inspect the distribution of the response variable</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb6"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb6-1"><a href="#cb6-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(d<span class="sc">$</span>Tmax, <span class="at">freq =</span> <span class="cn">FALSE</span>, <span class="at">breaks =</span> <span class="st">"Scott"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<p><img src="figures/hist_Tmax-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>The histogram suggests that the data is slightly left-skewed, with longer tails for temperatures below zero. This indicates that the commonly used normal distribution may not be the most appropriate choice for modeling daily maximum temperatures.</p>
<p>To address this, the <code>gamlss2</code> package provides the <a href="https://gamlss-dev.github.io/gamlss2/man/find_family.html"><code>find_family()</code></a> function, which helps identify the most suitable distribution by minimizing an information criterion, AIC by default. Here, we evaluate several continuous distributions available in the <a href="https://cran.r-project.org/package=gamlss.dist"><code>gamlss.dist</code></a> package.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb7"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb7-1"><a href="#cb7-1" aria-hidden="true" tabindex="-1"></a>fams <span class="ot">&lt;-</span> <span class="fu">find_family</span>(d<span class="sc">$</span>Tmax,</span>
<span id="cb7-2"><a href="#cb7-2" aria-hidden="true" tabindex="-1"></a>  <span class="at">families =</span> <span class="fu">c</span>(NO, TF, JSU, SEP4))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>.. NO family
.. .. IC = 92044.1 
.. TF family
.. .. IC = 92046.12 
.. JSU family
.. .. IC = 91975.24 
.. SEP4 family
.. .. IC = 91879.66 </code></pre>
</div>
<div class="sourceCode cell-code" id="cb9"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb9-1"><a href="#cb9-1" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(fams)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>      TF       NO      JSU     SEP4 
92046.12 92044.10 91975.24 91879.66 </code></pre>
</div>
</div>
<p>Here, the <code>SEP4</code> family appears to provide the best fit. To further assess its suitability, we can visualize the fitted density using</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb11"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb11-1"><a href="#cb11-1" aria-hidden="true" tabindex="-1"></a><span class="fu">fit_family</span>(d<span class="sc">$</span>Tmax, <span class="at">family =</span> SEP4)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>GAMLSS-RS iteration  1: Global Deviance = 91950.4801 eps = 0.044494     
GAMLSS-RS iteration  2: Global Deviance = 91907.9542 eps = 0.000462     
GAMLSS-RS iteration  3: Global Deviance = 91889.9353 eps = 0.000196     
GAMLSS-RS iteration  4: Global Deviance = 91880.9925 eps = 0.000097     
GAMLSS-RS iteration  5: Global Deviance = 91876.7482 eps = 0.000046     
GAMLSS-RS iteration  6: Global Deviance = 91874.399 eps = 0.000025     
GAMLSS-RS iteration  7: Global Deviance = 91872.4447 eps = 0.000021     
GAMLSS-RS iteration  8: Global Deviance = 91871.6641 eps = 0.000008     </code></pre>
</div>
<div class="cell-output-display">
<p><img src="figures/fit_family-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>After identifying a suitable distributional model, we can now incorporate covariates to estimate a full GAMLSS. Since temperature data exhibits a strong seasonal pattern, as illustrated in the following scatterplot</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb13"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb13-1"><a href="#cb13-1" aria-hidden="true" tabindex="-1"></a>d<span class="sc">$</span>yday <span class="ot">&lt;-</span> <span class="fu">as.POSIXlt</span>(d<span class="sc">$</span>date)<span class="sc">$</span>yday</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<p><img src="figures/season_scatter-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>It is essential to include a model term that captures these seasonal effects. The <code>gamlss2</code> package supports all model terms from the <a href="https://cran.r-project.org/package=mgcv"><code>mgcv</code></a> package, allowing us to use the <code>s()</code> constructor to model seasonality.</p>
<p>Additionally, we include a time trend to examine whether maximum temperatures have increased over the observed period. In the full GAMLSS model, each parameter of the selected <code>SEP4</code> distribution is estimated separately. To incorporate the time trend, we first create a new covariate, <code>year</code>, representing the long-term temporal effect</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb14"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb14-1"><a href="#cb14-1" aria-hidden="true" tabindex="-1"></a>d<span class="sc">$</span>year <span class="ot">&lt;-</span> <span class="fu">as.POSIXlt</span>(d<span class="sc">$</span>date)<span class="sc">$</span>year <span class="sc">+</span> <span class="dv">1900</span></span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>Next, we define the model formula for the four parameters of the <code>SEP4</code> distribution.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb15"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb15-1"><a href="#cb15-1" aria-hidden="true" tabindex="-1"></a>f <span class="ot">&lt;-</span> Tmax <span class="sc">~</span> <span class="fu">s</span>(yday, <span class="at">bs =</span> <span class="st">"cc"</span>, <span class="at">k =</span> <span class="dv">20</span>) <span class="sc">+</span> <span class="fu">s</span>(year) <span class="sc">+</span> <span class="fu">te</span>(yday, year, <span class="at">bs =</span> <span class="fu">c</span>(<span class="st">"cc"</span>, <span class="st">"cr"</span>)) <span class="sc">|</span></span>
<span id="cb15-2"><a href="#cb15-2" aria-hidden="true" tabindex="-1"></a>  <span class="fu">s</span>(yday, <span class="at">bs =</span> <span class="st">"cc"</span>, <span class="at">k =</span> <span class="dv">20</span>) <span class="sc">+</span> <span class="fu">s</span>(year) <span class="sc">+</span> <span class="fu">te</span>(yday, year, <span class="at">bs =</span> <span class="fu">c</span>(<span class="st">"cc"</span>, <span class="st">"cr"</span>)) <span class="sc">|</span></span>
<span id="cb15-3"><a href="#cb15-3" aria-hidden="true" tabindex="-1"></a>  <span class="fu">s</span>(yday, <span class="at">bs =</span> <span class="st">"cc"</span>, <span class="at">k =</span> <span class="dv">20</span>) <span class="sc">+</span> <span class="fu">s</span>(year) <span class="sc">+</span> <span class="fu">te</span>(yday, year, <span class="at">bs =</span> <span class="fu">c</span>(<span class="st">"cc"</span>, <span class="st">"cr"</span>)) <span class="sc">|</span></span>
<span id="cb15-4"><a href="#cb15-4" aria-hidden="true" tabindex="-1"></a>  <span class="fu">s</span>(yday, <span class="at">bs =</span> <span class="st">"cc"</span>, <span class="at">k =</span> <span class="dv">20</span>) <span class="sc">+</span> <span class="fu">s</span>(year) <span class="sc">+</span> <span class="fu">te</span>(yday, year, <span class="at">bs =</span> <span class="fu">c</span>(<span class="st">"cc"</span>, <span class="st">"cr"</span>))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>In this formula, the vertical bars <code>|</code> separate the specifications for each parameter of the <code>SEP4</code> distribution. The argument <code>bs = "cc"</code> specifies a cyclical spline to account for the seasonal effect, ensuring continuity at the beginning and end of the year, and argument <code>k</code> controls the dimension of the basis used to represent the smooth term.</p>
<p>Finally, we estimate the model using</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb16"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb16-1"><a href="#cb16-1" aria-hidden="true" tabindex="-1"></a>b <span class="ot">&lt;-</span> <span class="fu">gamlss2</span>(f, <span class="at">data =</span> d, <span class="at">family =</span> SEP4)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>GAMLSS-RS iteration  1: Global Deviance = 79918.3054 eps = 0.169526     
GAMLSS-RS iteration  2: Global Deviance = 79804.3251 eps = 0.001426     
GAMLSS-RS iteration  3: Global Deviance = 79759.9851 eps = 0.000555     
GAMLSS-RS iteration  4: Global Deviance = 79745.0554 eps = 0.000187     
GAMLSS-RS iteration  5: Global Deviance = 79738.3587 eps = 0.000083     
GAMLSS-RS iteration  6: Global Deviance = 79736.4781 eps = 0.000023     
GAMLSS-RS iteration  7: Global Deviance = 79735.2054 eps = 0.000015     
GAMLSS-RS iteration  8: Global Deviance = 79734.6366 eps = 0.000007     </code></pre>
</div>
</div>
<p>This approach allows us to flexibly capture both seasonal patterns and long-term trends in daily maximum temperatures.</p>
<p>After estimating the model, we can examine the model summary using</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb18"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb18-1"><a href="#cb18-1" aria-hidden="true" tabindex="-1"></a><span class="fu">summary</span>(b)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>Call:
gamlss2(formula = f, data = d, family = SEP4)
---
Family: SEP4 
Link function: mu = identity, sigma = log, nu = log, tau = log
*--------
Parameter: mu 
---
Coefficients:
            Estimate Std. Error t value Pr(&gt;|t|)    
(Intercept)  -0.9380     0.1507  -6.224    5e-10 ***
---
Smooth terms:
    s(yday) s(year) te(yday,year)
edf 16.3620  7.6634        17.997
*--------
Parameter: sigma 
---
Coefficients:
            Estimate Std. Error t value Pr(&gt;|t|)    
(Intercept)  1.88151    0.00144    1307   &lt;2e-16 ***
---
Smooth terms:
    s(yday) s(year) te(yday,year)
edf 15.9695  4.6506        5.9388
*--------
Parameter: nu 
---
Coefficients:
            Estimate Std. Error t value Pr(&gt;|t|)    
(Intercept)  0.63884    0.02157   29.62   &lt;2e-16 ***
---
Smooth terms:
    s(yday) s(year) te(yday,year)
edf 12.2077  6.0049        2.0948
*--------
Parameter: tau 
---
Coefficients:
            Estimate Std. Error t value Pr(&gt;|t|)    
(Intercept) 1.027653   0.009979     103   &lt;2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
---
Smooth terms:
    s(yday) s(year) te(yday,year)
edf  11.060   1.013        1.4077
*--------
n = 13665 df =  106.37 res.df =  13558.63
Deviance = 79734.6366 Null Dev. Red. = 13.21%
AIC = 79947.3743 elapsed = 20.67sec</code></pre>
</div>
</div>
<p>The summary output is structured similarly to those of <code>lm()</code> and <code>glm()</code>, with the key difference being that it provides results for all parameters of the selected distribution. Specifically, it displays the estimated linear coefficients (in this case, primarily the intercepts), along with the effective degrees of freedom for each smooth term. Additionally, the AIC and deviance values are reported.</p>
<p>To extract the AIC separately, we use</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb20"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb20-1"><a href="#cb20-1" aria-hidden="true" tabindex="-1"></a><span class="fu">AIC</span>(b)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>[1] 79947.37</code></pre>
</div>
</div>
<p>Similarly, the log-likelihood can be obtained with</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb22"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb22-1"><a href="#cb22-1" aria-hidden="true" tabindex="-1"></a><span class="fu">logLik</span>(b)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>'log Lik.' -39867.32 (df=106.3689)</code></pre>
</div>
<div class="sourceCode cell-code" id="cb24"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb24-1"><a href="#cb24-1" aria-hidden="true" tabindex="-1"></a><span class="fu">logLik</span>(b, <span class="at">newdata =</span> d)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>'log Lik.' -39867.32 (df=106.3689)</code></pre>
</div>
</div>
<p>Here we use the <code>newdata</code> argument just to show, that the log-likelihood can also be evaluated on, e.g., out-of-sample data.</p>
<p>Additionally, the estimated effects can be visualized instantly using</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb26"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb26-1"><a href="#cb26-1" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(b, <span class="at">which =</span> <span class="st">"effects"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<p><img src="figures/effects-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>This plot provides a direct visualization of the smooth effects included in the model, helping to interpret seasonal variations and long-term trends efficiently.</p>
<p>To assess the calibration of the estimated model, we examine the quantile residuals using a histogram, Q-Q plot, and worm plot.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb27"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb27-1"><a href="#cb27-1" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(b, <span class="at">which =</span> <span class="st">"resid"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<p><img src="figures/resids-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>These diagnostic plots indicate that the model is well-calibrated when using the <code>SEP4</code> distribution, demonstrating a good fit to the observed data.</p>
<p>Model predictions can be obtained for various statistical quantities, including the mean, quantiles, probability density function (PDF), and cumulative distribution function (CDF). To illustrate this, we first examine the marginal effect of the long-term time trend. For this purpose, we create a new data frame containing only the years of interest.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb28"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb28-1"><a href="#cb28-1" aria-hidden="true" tabindex="-1"></a>nd <span class="ot">&lt;-</span> <span class="fu">data.frame</span>(<span class="st">"year"</span> <span class="ot">=</span> <span class="dv">1981</span><span class="sc">:</span><span class="dv">2018</span>, <span class="st">"yday"</span> <span class="ot">=</span> <span class="dv">182</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>Next, we predict quantiles by first computing the estimated parameters</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb29"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb29-1"><a href="#cb29-1" aria-hidden="true" tabindex="-1"></a>par <span class="ot">&lt;-</span> <span class="fu">predict</span>(b, <span class="at">newdata =</span> nd)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>To compute, e.g., the 50% quantile (median), we extract the <a href="https://gamlss-dev.github.io/gamlss2/man/gamlss2.family.html"><code>gamlss2.family</code></a> of the fitted model and call the corresponding <code>$q()</code> (quantile) function provided by the family.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb30"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb30-1"><a href="#cb30-1" aria-hidden="true" tabindex="-1"></a>q50 <span class="ot">&lt;-</span> <span class="fu">family</span>(b)<span class="sc">$</span><span class="fu">q</span>(<span class="fl">0.5</span>, par)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>Similarly, we can compute the 10% and 90% quantiles</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb31"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb31-1"><a href="#cb31-1" aria-hidden="true" tabindex="-1"></a>q10 <span class="ot">&lt;-</span> <span class="fu">family</span>(b)<span class="sc">$</span><span class="fu">q</span>(<span class="fl">0.1</span>, par)</span>
<span id="cb31-2"><a href="#cb31-2" aria-hidden="true" tabindex="-1"></a>q90 <span class="ot">&lt;-</span> <span class="fu">family</span>(b)<span class="sc">$</span><span class="fu">q</span>(<span class="fl">0.9</span>, par)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<p>Finally, we visualize the long-term trend in temperature.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb32"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb32-1"><a href="#cb32-1" aria-hidden="true" tabindex="-1"></a><span class="fu">matplot</span>(nd<span class="sc">$</span>year, <span class="fu">cbind</span>(q10, q50, q90), <span class="at">type =</span> <span class="st">"l"</span>,</span>
<span id="cb32-2"><a href="#cb32-2" aria-hidden="true" tabindex="-1"></a>  <span class="at">lwd =</span> <span class="dv">2</span>, <span class="at">xlab =</span> <span class="st">"Year"</span>, <span class="at">ylab =</span> <span class="st">"Estimated Quantiles"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<p><img src="figures/long_term-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>The plot reveals an upward trend in the median temperature over time, highlighting the effects of long-term climate change.</p>
<p>To visualize exceedance probabilities for the 2019 season, we use the <code>$p()</code> function of the family object. For example, we can compute the probabilities of maximum temperatures exceeding 10, 11, 12, 13, and 14 °C as follows</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb33"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb33-1"><a href="#cb33-1" aria-hidden="true" tabindex="-1"></a>nd <span class="ot">&lt;-</span> <span class="fu">data.frame</span>(<span class="st">"year"</span> <span class="ot">=</span> <span class="dv">2019</span>, <span class="st">"yday"</span> <span class="ot">=</span> <span class="dv">0</span><span class="sc">:</span><span class="dv">365</span>)</span>
<span id="cb33-2"><a href="#cb33-2" aria-hidden="true" tabindex="-1"></a>par <span class="ot">&lt;-</span> <span class="fu">predict</span>(b, <span class="at">newdata =</span> nd)</span>
<span id="cb33-3"><a href="#cb33-3" aria-hidden="true" tabindex="-1"></a>Tmax <span class="ot">&lt;-</span> <span class="fu">rev</span>(<span class="fu">seq</span>(<span class="dv">0</span>, <span class="dv">14</span>, <span class="at">by =</span> <span class="dv">2</span>))</span>
<span id="cb33-4"><a href="#cb33-4" aria-hidden="true" tabindex="-1"></a>probs <span class="ot">&lt;-</span> <span class="fu">sapply</span>(Tmax, <span class="cf">function</span>(t) <span class="dv">1</span> <span class="sc">-</span> <span class="fu">family</span>(b)<span class="sc">$</span><span class="fu">p</span>(t, par))</span>
<span id="cb33-5"><a href="#cb33-5" aria-hidden="true" tabindex="-1"></a><span class="fu">colnames</span>(probs) <span class="ot">&lt;-</span> <span class="fu">paste0</span>(<span class="st">"Prob(Tmax &gt; "</span>, Tmax, <span class="st">")"</span>)</span>
<span id="cb33-6"><a href="#cb33-6" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(probs)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>     Prob(Tmax &gt; 14) Prob(Tmax &gt; 12) Prob(Tmax &gt; 10) Prob(Tmax &gt; 8)
[1,]    1.590061e-12    5.501117e-10    7.673712e-08   4.605888e-06
[2,]    1.635692e-12    5.493598e-10    7.520458e-08   4.469749e-06
[3,]    1.630585e-12    5.369963e-10    7.269465e-08   4.302274e-06
[4,]    1.585621e-12    5.161699e-10    6.952153e-08   4.115322e-06
[5,]    1.514233e-12    4.901235e-10    6.598375e-08   3.919732e-06
[6,]    1.428635e-12    4.618508e-10    6.234650e-08   3.725021e-06
     Prob(Tmax &gt; 6) Prob(Tmax &gt; 4) Prob(Tmax &gt; 2) Prob(Tmax &gt; 0)
[1,]   0.0001275031    0.001753963     0.01299744     0.05663021
[2,]   0.0001233908    0.001701366     0.01267788     0.05562356
[3,]   0.0001188968    0.001647312     0.01236082     0.05464662
[4,]   0.0001142025    0.001592977     0.01204944     0.05370176
[5,]   0.0001094721    0.001539443     0.01174679     0.05279154
[6,]   0.0001048512    0.001487695     0.01145582     0.05191868</code></pre>
</div>
</div>
<p>To illustrate these exceedance probabilities, we plot them over the course of the year</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb35"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb35-1"><a href="#cb35-1" aria-hidden="true" tabindex="-1"></a>col <span class="ot">&lt;-</span> colorspace<span class="sc">::</span><span class="fu">heat_hcl</span>(<span class="fu">ncol</span>(probs))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<p><img src="figures/vis_probs-1.png" class="img-fluid" width="672"></p>
</div>
</div>
<p>The plot reveals that even at this high-altitude station, the probability of <code>Tmax</code> &gt; 14°C reaches approximately 5% during summer. This is particularly striking considering that the Zugspitze once had a permanent glacier field, emphasizing the impact of rising temperatures in this region. Likewise, the probability of <code>Tmax</code> &gt; 0°C during the winter months is also about 5%, highlighting significant temperature patterns.</p>
</section>
</section>

</main>
<!-- /main column -->
<script id="quarto-html-after-body" type="application/javascript">
window.document.addEventListener("DOMContentLoaded", function (event) {
  const toggleBodyColorMode = (bsSheetEl) => {
    const mode = bsSheetEl.getAttribute("data-mode");
    const bodyEl = window.document.querySelector("body");
    if (mode === "dark") {
      bodyEl.classList.add("quarto-dark");
      bodyEl.classList.remove("quarto-light");
    } else {
      bodyEl.classList.add("quarto-light");
      bodyEl.classList.remove("quarto-dark");
    }
  }
  const toggleBodyColorPrimary = () => {
    const bsSheetEl = window.document.querySelector("link#quarto-bootstrap");
    if (bsSheetEl) {
      toggleBodyColorMode(bsSheetEl);
    }
  }
  toggleBodyColorPrimary();  
  const icon = "";
  const anchorJS = new window.AnchorJS();
  anchorJS.options = {
    placement: 'right',
    icon: icon
  };
  anchorJS.add('.anchored');
  const isCodeAnnotation = (el) => {
    for (const clz of el.classList) {
      if (clz.startsWith('code-annotation-')) {                     
        return true;
      }
    }
    return false;
  }
  const clipboard = new window.ClipboardJS('.code-copy-button', {
    text: function(trigger) {
      const codeEl = trigger.previousElementSibling.cloneNode(true);
      for (const childEl of codeEl.children) {
        if (isCodeAnnotation(childEl)) {
          childEl.remove();
        }
      }
      return codeEl.innerText;
    }
  });
  clipboard.on('success', function(e) {
    // button target
    const button = e.trigger;
    // don't keep focus
    button.blur();
    // flash "checked"
    button.classList.add('code-copy-button-checked');
    var currentTitle = button.getAttribute("title");
    button.setAttribute("title", "Copied!");
    let tooltip;
    if (window.bootstrap) {
      button.setAttribute("data-bs-toggle", "tooltip");
      button.setAttribute("data-bs-placement", "left");
      button.setAttribute("data-bs-title", "Copied!");
      tooltip = new bootstrap.Tooltip(button, 
        { trigger: "manual", 
          customClass: "code-copy-button-tooltip",
          offset: [0, -8]});
      tooltip.show();    
    }
    setTimeout(function() {
      if (tooltip) {
        tooltip.hide();
        button.removeAttribute("data-bs-title");
        button.removeAttribute("data-bs-toggle");
        button.removeAttribute("data-bs-placement");
      }
      button.setAttribute("title", currentTitle);
      button.classList.remove('code-copy-button-checked');
    }, 1000);
    // clear code selection
    e.clearSelection();
  });
  function tippyHover(el, contentFn) {
    const config = {
      allowHTML: true,
      content: contentFn,
      maxWidth: 500,
      delay: 100,
      arrow: false,
      appendTo: function(el) {
          return el.parentElement;
      },
      interactive: true,
      interactiveBorder: 10,
      theme: 'quarto',
      placement: 'bottom-start'
    };
    window.tippy(el, config); 
  }
  const noterefs = window.document.querySelectorAll('a[role="doc-noteref"]');
  for (var i=0; i<noterefs.length; i++) {
    const ref = noterefs[i];
    tippyHover(ref, function() {
      // use id or data attribute instead here
      let href = ref.getAttribute('data-footnote-href') || ref.getAttribute('href');
      try { href = new URL(href).hash; } catch {}
      const id = href.replace(/^#\/?/, "");
      const note = window.document.getElementById(id);
      return note.innerHTML;
    });
  }
      let selectedAnnoteEl;
      const selectorForAnnotation = ( cell, annotation) => {
        let cellAttr = 'data-code-cell="' + cell + '"';
        let lineAttr = 'data-code-annotation="' +  annotation + '"';
        const selector = 'span[' + cellAttr + '][' + lineAttr + ']';
        return selector;
      }
      const selectCodeLines = (annoteEl) => {
        const doc = window.document;
        const targetCell = annoteEl.getAttribute("data-target-cell");
        const targetAnnotation = annoteEl.getAttribute("data-target-annotation");
        const annoteSpan = window.document.querySelector(selectorForAnnotation(targetCell, targetAnnotation));
        const lines = annoteSpan.getAttribute("data-code-lines").split(",");
        const lineIds = lines.map((line) => {
          return targetCell + "-" + line;
        })
        let top = null;
        let height = null;
        let parent = null;
        if (lineIds.length > 0) {
            //compute the position of the single el (top and bottom and make a div)
            const el = window.document.getElementById(lineIds[0]);
            top = el.offsetTop;
            height = el.offsetHeight;
            parent = el.parentElement.parentElement;
          if (lineIds.length > 1) {
            const lastEl = window.document.getElementById(lineIds[lineIds.length - 1]);
            const bottom = lastEl.offsetTop + lastEl.offsetHeight;
            height = bottom - top;
          }
          if (top !== null && height !== null && parent !== null) {
            // cook up a div (if necessary) and position it 
            let div = window.document.getElementById("code-annotation-line-highlight");
            if (div === null) {
              div = window.document.createElement("div");
              div.setAttribute("id", "code-annotation-line-highlight");
              div.style.position = 'absolute';
              parent.appendChild(div);
            }
            div.style.top = top - 2 + "px";
            div.style.height = height + 4 + "px";
            let gutterDiv = window.document.getElementById("code-annotation-line-highlight-gutter");
            if (gutterDiv === null) {
              gutterDiv = window.document.createElement("div");
              gutterDiv.setAttribute("id", "code-annotation-line-highlight-gutter");
              gutterDiv.style.position = 'absolute';
              const codeCell = window.document.getElementById(targetCell);
              const gutter = codeCell.querySelector('.code-annotation-gutter');
              gutter.appendChild(gutterDiv);
            }
            gutterDiv.style.top = top - 2 + "px";
            gutterDiv.style.height = height + 4 + "px";
          }
          selectedAnnoteEl = annoteEl;
        }
      };
      const unselectCodeLines = () => {
        const elementsIds = ["code-annotation-line-highlight", "code-annotation-line-highlight-gutter"];
        elementsIds.forEach((elId) => {
          const div = window.document.getElementById(elId);
          if (div) {
            div.remove();
          }
        });
        selectedAnnoteEl = undefined;
      };
      // Attach click handler to the DT
      const annoteDls = window.document.querySelectorAll('dt[data-target-cell]');
      for (const annoteDlNode of annoteDls) {
        annoteDlNode.addEventListener('click', (event) => {
          const clickedEl = event.target;
          if (clickedEl !== selectedAnnoteEl) {
            unselectCodeLines();
            const activeEl = window.document.querySelector('dt[data-target-cell].code-annotation-active');
            if (activeEl) {
              activeEl.classList.remove('code-annotation-active');
            }
            selectCodeLines(clickedEl);
            clickedEl.classList.add('code-annotation-active');
          } else {
            // Unselect the line
            unselectCodeLines();
            clickedEl.classList.remove('code-annotation-active');
          }
        });
      }
  const findCites = (el) => {
    const parentEl = el.parentElement;
    if (parentEl) {
      const cites = parentEl.dataset.cites;
      if (cites) {
        return {
          el,
          cites: cites.split(' ')
        };
      } else {
        return findCites(el.parentElement)
      }
    } else {
      return undefined;
    }
  };
  var bibliorefs = window.document.querySelectorAll('a[role="doc-biblioref"]');
  for (var i=0; i<bibliorefs.length; i++) {
    const ref = bibliorefs[i];
    const citeInfo = findCites(ref);
    if (citeInfo) {
      tippyHover(citeInfo.el, function() {
        var popup = window.document.createElement('div');
        citeInfo.cites.forEach(function(cite) {
          var citeDiv = window.document.createElement('div');
          citeDiv.classList.add('hanging-indent');
          citeDiv.classList.add('csl-entry');
          var biblioDiv = window.document.getElementById('ref-' + cite);
          if (biblioDiv) {
            citeDiv.innerHTML = biblioDiv.innerHTML;
          }
          popup.appendChild(citeDiv);
        });
        return popup.innerHTML;
      });
    }
  }
});
</script>
</div> <!-- /content -->



</body></html>