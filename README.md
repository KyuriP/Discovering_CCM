# MSBBSS Thesis Project: "Discovering Cyclic Causal Models in Psychological Research"

This repository contains all the files concerning the thesis project.

Checklist
================

Content Overview
- <a href="#introduction" id="toc-introduction"><span
  class="toc-section-number">1</span> Introduction</a>
  - <a href="#background" id="toc-background"><span
    class="toc-section-number">1.1</span> Background</a>
    - <a href="#graphical-models" id="toc-graphical-models"><span
      class="toc-section-number">1.1.1</span> Graphical Models</a> :heavy_check_mark:
    - <a href="#acylic-vs.-cyclic-graphs" 
      id="toc-acylic-vs.-cyclic-graphs"><span
      class="toc-section-number">1.1.2</span> Acylic vs. Cyclic Graphs:</a> :heavy_check_mark:
    - <a href="#primer-on-constraint-based-methods"
      id="toc-primer-on-constraint-based-methods"><span
      class="toc-section-number">1.1.3</span> Primer on Constraint-based
      Methods</a> :heavy_check_mark:
- <a href="#methods" id="toc-methods"><span
  class="toc-section-number">2</span> Methods</a>
  - <a href="#ccd" id="toc-ccd"><span class="toc-section-number">2.1</span>
    CCD</a> :heavy_check_mark:
  - <a href="#fci" id="toc-fci"><span class="toc-section-number">2.2</span>
    FCI</a> :heavy_check_mark:
  - <a href="#cci" id="toc-cci"><span class="toc-section-number">2.3</span>
    CCI</a> :heavy_check_mark:
  - <a href="#simulation" id="toc-simulation"><span
    class="toc-section-number">2.4</span> Simulation</a> :heavy_check_mark:
  - <a href="#emprical-example" id="toc-emprical-example"><span
    class="toc-section-number">2.5</span> Emprical Example</a>
- <a href="#results" id="toc-results"><span
  class="toc-section-number">3</span> Results</a>
  - <a href="#simulation-study" id="toc-simulation-study"><span
    class="toc-section-number">3.1</span> Simulation study
    <input type="checkbox" name="simulation"/></a> 
  - <a href="#empirical-application" id="toc-empirical-application"><span
    class="toc-section-number">3.2</span> Empirical application
    <input type="checkbox" name="empirical"/></a>
- <a href="#discussion" id="toc-discussion"><span
  class="toc-section-number">4</span> Discussion</a>
  - <a href="#limitations" id="toc-limitations"><span
    class="toc-section-number">4.1</span> Limitations
    <input type="checkbox" name="limitations"/></a>
  - <a href="#future-research" id="toc-future-research"><span
    class="toc-section-number">4.2</span> Future research
    <input type="checkbox" name="future"/></a>

<hr>

# To-dos

<table style="width:99%;">
<caption>Tentative Schedule</caption>
<colgroup>
<col style="width: 19%" />
<col style="width: 25%" />
<col style="width: 53%" />
</colgroup>
<thead>
<tr class="header">
<th>WEEK</th>
<th>TO-DO</th>
<th>DETAIL</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td> <s>week4</s> </td>
<td>Work on simulation</td>
<td><ul>
<li>Run CCD, FCI, CCI on an example model.</li>
<li>Update the code.</li>
</ul></td>
</tr>
<tr class="even">
<td><s>week5</s></td>
<td>Work on simulation</td>
<td><ul>
<li>Run CCD, FCI, CCI on an example model.</li>
<li>Update the code.</li>
</ul></td>
</tr>
<tr class="odd">
<td><s>week6</s></td>
<td><p>Work on simulation</p>
<p>Research seminar</p></td>
<td><ul>
<li>Trouble-shooting code.</li>
<li>Create summary figures / tables.</li>
<li>Seminar about creating an archive.</li>
</ul></td>
</tr>
<tr class="even">
<td><s>week7</s></td>
<td>Work on simulation &amp; Consultancy session</td>
<td><ul>
<li>Trouble-shooting code.</li>
<li>Create summary figures / tables.</li>
</ul></td>
</tr>
<tr class="odd">
<td>week8</td>
<td>Work on simulation &amp; Consultancy session</td>
<td><ul>
<li>Trouble-shooting code.</li>
<li>Create summary figures / tables.</li>
</ul></td>
</tr>
<tr class="even">
<td>week9</td>
<td>Work on simulation</td>
<td><ul>
<li>Review the results.</li>
<li>Update the simulation design and evaluation metrics.</li>
</ul></td>
</tr>
<tr class="odd">
<td>week10</td>
<td>Work on simulation &amp; update Methods section</td>
<td><ul>
<li>Update the results section.</li>
<li>Update the overview of algorithms and description of each
algorithm.</li>
</ul></td>
</tr>
<tr class="even">
<td>week11</td>
<td>Work on empirical analysis &amp; Consultancy session</td>
<td><ul>
<li>Choose an empirical data.</li>
<li>Try running the algorithm(s).</li>
</ul></td>
</tr>
<tr class="odd">
<td>week12</td>
<td>Work on empirical analysis &amp; Consultancy session & Peer feedback</td>
<td><ul>
<li>Try running the algorithm(s).</li>
<li>Trouble-shooting the code.</li>
<li>Send a draft (peer feedback).</li>
</ul></td>
</tr>
<tr class="even">
<td>week13</td>
<td>Work on empirical analysis &amp; Discussion</td>
<td><ul>
<li>Review the empirical results.</li>
<li>Update the results section.</li>
<li>Think of limitations/future extensions.</li>
</ul></td>
</tr>
<tr class="odd">
<td>week14</td>
<td>Update Results &amp; Discussion section</td>
<td><ul>
<li>Write a draft.</li>
<li>Ask for feedback.</li>
</ul></td>
</tr>
<tr class="even">
<td>week15</td>
<td>Update Methods, Results, and Discussion section</td>
<td><ul>
<li>Rewrite based on the feedback.</li>
</ul></td>
</tr>
<tr class="odd">
<td>week16</td>
<td>General Review</td>
<td><ul>
<li>Write the abstract.</li>
<li>Ask for feedback.</li>
</ul></td>
</tr>
<tr class="even">
<td>week17</td>
<td>General Review &amp; Research seminar</td>
<td><ul>
<li>Review the whole paper.</li>
<li>(Do what has not been done yet.)</li>
<li>Seminar about writing the abstract</li>
</ul></td>
</tr>
<tr class="odd">
<td>week18</td>
<td>Final check &amp; Work on Archive</td>
<td><ul>
<li>Final review.</li>
<li>Make sure the archive is more or less ready.</li>
</ul></td>
</tr>
<tr class="even">
<td>week19</td>
<td><p class="r">
Thesis Deadline
</p>
<p>&amp; Work on slides</p></td>
<td><ul>
<li>Submit the thesis (May 8th, 15:00).</li>
<li>Check the archive.</li>
<li>Make the slides for the defense.</li>
</ul></td>
</tr>
<tr class="odd">
<td>week20</td>
<td><p class="r">
Archive Deadline
</p>
<p>&amp; Work on poster</p></td>
<td><ul>
<li>Complete the archive (May 15th, 15:00).</li>
<li>Check the slides for defense.</li>
<li>Make the poster.</li>
</ul></td>
</tr>
<tr class="even">
<td>week21</td>
<td><p class="r">
Thesis Defense
</p></td>
<td><ul>
<li>Defense (May 22nd or 23rd)</li>
<li>Good luck.</li>
<li>Check the poster.</li>
</ul></td>
</tr>
<tr class="odd">
<td>week22</td>
<td>Poster Deadline</td>
<td><ul>
<li>Join the poster fair (May 30th)</li>
<li>The end.</li>
</ul></td>
</tr>
</tbody>
</table>
