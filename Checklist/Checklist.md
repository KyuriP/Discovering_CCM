---
title: "Checklist"
subtitle: 'Planning of the thesis project "Discivering Cyclic Causal Models in Psychological Research"'
author: "Kyuri Park"
date: "January 16, 2023"
format: 
  html:
    embed-resources: true
toc: true
toccolor: "#4169E1"
css: style.css
editor: visual
number-sections: true
keep-md: true
title-block-banner: true
title-block-banner-color: "#F5FFFA"
---


<hr>
# Content Overview {.unnumbered .sectitle}

# Introduction

-   Introduce the network theory <input type="checkbox" name="intronetwork"checked/>
-   Introduce causal discovery <input type="checkbox" name="intro-causaldiscovery"checked/>
-   Introduce the problem --\> What are we going to do? (RQ) <input type="checkbox" name="introprobs" checked/>
-   Connection to the network theory/modeling? <input type="checkbox"/>

## Background

### Graphical Models

-   Provide the basic concepts of graphical models <input type="checkbox" name="intro-graph" checked/>

### Acylic vs. Cyclic Graphs:

-   Explain the required assumptions for constraint-based methods <input type="checkbox" name="intro-assump" checked/>
-   Illustrate difficulties in estimating cyclic graphs compared to DAGs <input type="checkbox" name="cycleprob" checked/>

### Primer on Constraint-based Methods

-   Explain how constraint based methods work <input type="checkbox" name="primer" checked/>

# Methods

-   Introduce three algorithms we are going to examine <input type="checkbox" name="intro-algo" checked/>
-   Table of algorithm overview (assumptions / output representation) and explain the differences <input type="checkbox" name="overview-algo"/>
-   Introduce empirical data (McNally 2017?) <input type="checkbox" name="mcnally"/>

## CCD

-   Output representation: Markov equivalence class of DCGs (PAG) <input type="checkbox" name="ccdoutput" checked/>
-   CCD trace <input type="checkbox" name="ccdtrace" checked/>

## FCI

-   Output representation: Markov equivalence class of MAG (PAG) <input type="checkbox" name="fcioutput"/>
-   FCI trace <input type="checkbox" name="fcitrace"/>

## CCI

-   Output representation: Markov equivalence class of MAAG (PAAG?) <input type="checkbox" name="ccioutput"/>
-   CCI trace <input type="checkbox" name="ccitrace"/>

## Simulation

-   Simulation Design <input type="checkbox" name="simdesign" checked/>
-   Data Generating Process <input type="checkbox" name="datageneration" checked/>
-   Evaluation Metrics <input type="checkbox" name="evalmetrics" checked/>
-   Needs to describe them all in more detail <input type="checkbox"/>

## Emprical Example

-   Describe the example data in detail. <input type="checkbox" name="intro_empirical"/>

# Results

## Simulation study <input type="checkbox" name="simulation"/>

## Empirical application <input type="checkbox" name="empirical"/>

# Discussion

## Limitations <input type="checkbox" name="limitations"/>

## Future research <input type="checkbox" name="future"/>

<hr>

# To-dos {.unnumbered .sectitle}

+---------------+--------------------+------------------------------------------+
| WEEK          | TO-DO              | DETAIL                                   |
+===============+====================+==========================================+
| week4         | Work on simulation | - Run CCD, FCI, CCI on an example model. |
|               |                    | - Update the code.                       |
+---------------+--------------------+------------------------------------------+
| week5         | Work on simulation | - Run CCD, FCI, CCI on an example model. |
|               |                    | - Update the code.                       |
+---------------+--------------------+------------------------------------------+
| week6         | Work on simulation | - Trouble-shooting code.                 |
|               |                    | - Create summary figures / tables.       |
|               | Research seminar   | - Seminar about creating an archive.     |
+---------------+--------------------+------------------------------------------+
| week7         | Work on simulation | - Trouble-shooting code.                 |
|               | & Consultancy      | - Create summary figures / tables.       |
+---------------+--------------------+------------------------------------------+
| week8         | Work on simulation | - Trouble-shooting code.                 |
|               | & Consultancy      | - Create summary figures / tables.       |
+---------------+--------------------+------------------------------------------+
| week9         | Work on simulation | - Review the results.                    |
|               |                    | - Update the simulation design and       |
|               |                    |                       evaluation metrics.|
+---------------+--------------------+------------------------------------------+
| week10        | Work on simulation | - Update the results section.            |
|               | & update Methods   | - Update the overview of algorithms      |
|               |  section           |   and description of each algorithm.     |
+---------------+--------------------+------------------------------------------+
| week11        | Work on empirical  | - Choose an empirical data.              |
|               | analysis           | - Try running the algorithm(s).          |
|               | & Consultancy      |                                          |
+---------------+--------------------+------------------------------------------+
| week12        | Work on empirical  | - Try running the algorithm(s).          |
|               | analysis           | - Trouble-shooting the code.             |
|               | & Consultancy      | - Send a draft (peer feedback)           |
+---------------+--------------------+------------------------------------------+
| week13        | Work on empirical  | - Review the empirical results.          |
|               | analysis &         | - Update the results section.            |
|               | Discussion         | - Think of limitations/future extensions.|
+---------------+--------------------+------------------------------------------+
| week14        | Update Results &   | - Write a draft.                         |
|               | Discussion section | - Ask for feedback.                      |
+---------------+--------------------+------------------------------------------+
| week15        | Update Methods,    | - Rewrite based on the feedback.         |
|               | Results, and       |                                          |
|               | Discussion section |                                          |
+---------------+--------------------+------------------------------------------+
| week16        | General Review     | - Write the abstract.                    |
|               |                    | - Ask for feedback.                      |
+---------------+--------------------+------------------------------------------+
| week17        | General Review     | - Review the whole paper.                |
|               |          &         | - (Do what has not been done yet.)       |
|               | Research seminar   | - Seminar about writing the abstract     |
+---------------+--------------------+------------------------------------------+
| week18        | Final check &      | - Final review.                          |
|               | Work on Archive    | - Make sure the archive is more or less 
|               |                    |                                 ready.   |
+---------------+--------------------+------------------------------------------+
| week19        |<p class="r">Thesis | - Submit the thesis (May 8th, 15:00).    |
|               | Deadline</p>  &    | - Check the archive.                     |
|               |   Work on slides   | - Make the slides for the defense.       |
+---------------+--------------------+------------------------------------------+
| week20        |<p class="r">Archive|                                          |
|               |Deadline</p>        | - Complete the archive (May 15th, 15:00).|
|               | & Work on poster   | - Check the slides for defense.          |
|               |                    | - Make the poster.                       |
+---------------+--------------------+------------------------------------------+
| week21        |<p class="r">Thesis | - Defense (May 22nd or 23rd)             |       
|               | Defense</p>        | - Good luck.                             |
|               |                    | - Check the poster.                      |
+---------------+--------------------+------------------------------------------+
| week22        | Poster Deadline    | - Join the poster fair (May 30th)        |
|               |                    | - The end.                               |
+---------------+--------------------+------------------------------------------+

: Tentative Schedule