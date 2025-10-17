# *--------------------------------------------------------------------------* #
# * Configuration                                                            * #
# *--------------------------------------------------------------------------* #
from snakemake.utils import validate


include: "utils.smk"
include: "constants.smk"


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.json")


pepfile: "config/pep/config.yaml"


if config["dir_run"] and config["dir_run"] is not None:

    workdir: config["dir_run"]


# *--------------------------------------------------------------------------* #
# * Constant-like variables                                                  * #
# *--------------------------------------------------------------------------* #
SAMPLES = pep.sample_table["sample_name"]
SPECIES = config["species"]
CALLERS = sorted(config["callers"])
MAPPER = config["mapper"]
ANNOTATORS = config["annotators"]


# *--------------------------------------------------------------------------* #
# * Wildcard constraints                                                     * #
# *--------------------------------------------------------------------------* #
wildcard_constraints:
    sample=r"|".join(SAMPLES),
    type_sv=r"|".join(TYPES_SV),
    caller=r"|".join(CALLERS),


# *--------------------------------------------------------------------------* #
# * Files and directories required by rules                                  * #
# *--------------------------------------------------------------------------* #
path_cache_snpeff = (
    f"{config['cache_snpeff']}/{config['genome']}.{config['version_snpeff']}"
)
path_cache_vep = f"{config['cache_vep']}/{config['species']}/{config['version_vep']}_{config['genome']}"

vcfs_svision = multiext(
    "svision/{sample}/chroms/{sample}",
    *[f".{chrom}.svision.s{config['min_reads']}.graph.vcf" for chrom in CHROMS],
)
