Models
======

This folder contains the model used in the tutorials, as well as some
pre-curated thermodynamic information. The model that has such pre-curated
thermodynamic information will have a subfolder to its name, in which the
information is stored.

Pre-curated models
------------------

This model has been pre-curated for thermodynamics analysis, which means it already has the necessary annotations to be mapped to our thermodynamic data. Here is a short example for ``ipbe.mat``:


.. code-block:: matlab

    from cobra.io import load_matlab_model

    from pytfa.io import load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data


    # Load the model
    cobra_model = load_matlab_model('../models/small_ecoli.mat')
    
    # Load reaction DB
    thermo_data = load_thermoDB('../data/thermo_data.thermodb')
    lexicon = read_lexicon('../models/small_ecoli/lexicon.csv')
    compartment_data = read_compartment_data('../models/small_ecoli/compartment_data.json')

    # Initialize the model
    tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
    tmodel.name = 'tutorial'
    
    # Annotate the model
    annotate_from_lexicon(tmodel, lexicon)
    apply_compartment_data(tmodel, compartment_data)

    ## TFA conversion
    tmodel.prepare()
    tmodel.convert()

Publications
------------

The model ipbe is published together with the phenomapping repository in the publication:

+-----------+-------------------------------------------------------------------------------+
| ipbe.mat  | Orth, J. D., Conrad, T. M., Na, J., Lerman, J. A., Nam, H., Feist, A. M.,     |
|           | Palsson, & B. Ø. (2011). A comprehensive genome‐scale reconstruction of       |
|           | Escherichia coli metabolism—2011.    											|
|           | Molecular Systems Biology, 7(1):535, doi:10.1038/msb.2011.65.                 |
+-----------+-------------------------------------------------------------------------------+

  


