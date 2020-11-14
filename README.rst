PhenoMapping
============
|Documentation Status| |Build Status| |Codecov| |Codacy branch grade| |license|


**PhenoMapping**: A workflow for genome-scale models to

(1) analyse the context-specific metabolic function: at minimal media, with thermodynamic constraints, with metabolomics data integrated, with transcriptomics data integrated

(2) study essentiality at those conditions

(3) map essentiality to the underlying context-specific condition (active constraints)



This code is the release for the study of life-stage-specific metabolic function in the malaria parasite Plasmodium berghei.

It was applied to analyse high-throughput gene knockout data in the blood and liver stage development of the malaria parasite **Plasmodium berghei** using the genome-scale model of this organism iPbe (doi.org/10.1016/j.cell.2019.10.030).

It was also applied to analyse high-throughput gene knockout data for tachyzoites of **Toxoplasma gondii** using the genome-scale model of this organism iTgo (doi.org/10.1016/j.chom.2020.01.002).



You will need matTFA_ to run it and TEX-FBA_ to perform the transcriptomics related analyses.

We recommend the stable combination of **MATLAB** (any version between 2016a and 2019a) and **CPLEX** 12.7 (also freely downloadable from the IBM Academic initiative) to run PhenoMapping.


.. _Manuscript: Stanway R. R., Bushell E., Chiappino-Pepe A., Roques M., Sanderson T., Franke-Fayard B., Caldelari R., Golomingi M., Nyonda M., Pandey V., Schwach F., Chevalley S., Ramesar J., Metcalf T., Herd C., Burda P. C., Rayner J. C., Soldati-Favre D., Janse C., Hatzimanikatis V., Billker O., Heussler V. T, Cell (2019). Genome-Scale Identification of Essential Metabolic Processes for Targeting the Plasmodium Liver Stage. https://doi.org/10.1016/j.cell.2019.10.030.
.. _matTFA: https://github.com/EPFL-LCSB/matTFA
.. _TEX-FBA: https://github.com/EPFL-LCSB/texfba
.. _Documentation: https://phenomapping.readthedocs.io/en/latest/solver.html
.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/phenomapping/blob/master/LICENSE.txt
.. |Documentation Status| image:: https://readthedocs.org/projects/phenomapping/badge/?version=latest
   :target: http://phenomapping.readthedocs.io/en/latest/?badge=latest
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/phenomapping.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/phenomapping
.. |Codecov| image:: https://img.shields.io/codecov/c/github/EPFL-LCSB/phenomapping.svg
   :target: https://codecov.io/gh/EPFL-LCSB/phenomapping
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/46bab484396946a8be07a82276f3e9dc/master.svg
   :target: https://www.codacy.com/app/realLCSB/phenomapping


Quick-start
============

First clone or fork this repository at your desired location or path and navigate to the repository's top directory: 

.. code:: bash

    git clone https://github.com/EPFL-LCSB/phenomapping.git
    cd phenomapping

If not done before, install CPLEX to get started with PhenoMapping.

A fully running example of PhenoMapping is available at tutorials/tutorial_basics.mat


Description of files in this repository
---------------------------------------
**tutorials/tutorial_basics.m** - Start here for explanations and examples on preparing and loading the model, preparing data, and running modules.

**tests** - Contains all Matlab scripts to run independently the PhenoMapping modules.

**tests/settings.m** - Template to adapt to your model. See adapted template for iPbe (settings_ipbeliver.m, settings_ipbeblood.m) and iTgo (settings_itgo.m)

**tests/test_core_modulename.m** - Script of each module. It is recommended to run these individually and separately in the order defined in tutorial_basics.

**tests/ref/pbe** - Contains .mat files with the data integrated into iPbe in Stanway et al. and used for the example case in tutorials/tutorial_basics.mat

**models/** - Contains the non-context specific model iPbe used in Stanway et al. and used for the example case in tutorials/tutorial_basics.mat

**phenomapping/** - (subfolder) Contains all functions required to run PhenoMapping



Please, let us know if you have any suggestion, comment, or problem. We will be happy to discuss and help.

Contact: Dr. Anush Chiappino-Pepe (anush.chiappinopepe@alumni.epfl.ch)



License
=======
The code in this repository is licensed under the terms of APACHEv2.0 as specified by the LICENSE `<https://github.com/EPFL-LCSB/phenomapping/blob/master/LICENSE.txt>`_ file



References
==========
Please cite the following reference for the PhenoMapping package:

Stanway R. R., Bushell E., Chiappino-Pepe A., Roques M., Sanderson T., Franke-Fayard B., Caldelari R., Golomingi M., Nyonda M., Pandey V., Schwach F., Chevalley S., Ramesar J., Metcalf T., Herd C., Burda P. C., Rayner J. C., Soldati-Favre D., Janse C., Hatzimanikatis V., Billker O., Heussler V. T, Cell (2019). Genome-Scale Identification of Essential Metabolic Processes for Targeting the Plasmodium Liver Stage. https://doi.org/10.1016/j.cell.2019.10.030.



Also used in
============

Krishnan A., Kloehn J., Lunghi M., Chiappino-Pepe A., Waldman B.S., Nicolas D., Varesio E., Hehl A., Lourido S., Hatzimanikatis V., Soldati-Favre S., (2020). Functional and Computational Genomics Reveal Unprecedented Flexibility in Stage-Specific Toxoplasma Metabolism. Cell Host & Microbe 27(2), 290-306.e11. https://doi.org/10.1016/j.chom.2020.01.002

