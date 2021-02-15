************
PanACoTA API
************

Hereafter is the PanACoTA API, describing all modules and their functions.

``PanACoTA`` package contains:

    - 2 :doc:`submodules <PanACoTA.utils>`: ``utils`` and ``utils_pangenome``
    - 1 subpackage, called :doc:`subcommands <PanACoTA.subcommands>`, containing all subcommands main scripts
    - 6 subpackages, corresponding to the 6 subcommands:
        * :doc:`prepare_module<PanACoTA.prepare_module>`
        * :doc:`annotate_module<PanACoTA.annotate_module>`
        * :doc:`pangenome_module<PanACoTA.pangenome_module>`,
        * :doc:`corepers_module<PanACoTA.corepers_module>`,
        * :doc:`align_module<PanACoTA.align_module>`,
        * :doc:`tree_module<PanACoTA.tree_module>`.
    - the module to run all subcommands:
        :doc:`all_module<PanACoTA.all_modules>`.


.. toctree::
    :hidden:

    PanACoTA.utils
    PanACoTA.subcommands
    PanACoTA.prepare_module
    PanACoTA.annotate_module
    PanACoTA.pangenome_module
    PanACoTA.corepers_module
    PanACoTA.align_module
    PanACoTA.tree_module