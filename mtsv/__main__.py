import click
import yaml
import subprocess as sp
from collections import OrderedDict
from mtsv import VERSION, PROG_NAME, CONFIG_STRING, CLUSTER_CONFIG
from mtsv.utils import (
    snake_path, error, warn, info, data_path, specfile_path)
from mtsv.scripts.datastore_management import (
    get_datastore_info_table, remove_datastore)

SPECFILE_PATH = specfile_path("validation")

SNAKEPATH = snake_path("Snakefile")

# Args that change behavior of snakemake
PASSTHROUGH_ARGS = [
    "-S", "--summary", "--cleanup-shadow", "--delete-all-output",
    "--delete-temp-output", "--report", "--export-cwl",
    "--list", "-l", "--list-target-rules", "--lt", "--dag", "--rulegraph",
    "--filegraph", "--d3dag", "--detailed-summary", "-D", "--archive",
    "--cleanup-metadata"]


class OrderedGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super(OrderedGroup, self).__init__(name, commands, **attrs)
        #: the registered subcommands by their exported names.
        self.commands = commands or OrderedDict()

    def list_commands(self, ctx):
        return self.commands

def run_subprocess(cmd):
    try:
        result = sp.run(cmd, stdout=sp.PIPE, check=True)
    except sp.CalledProcessError as e:
        error(
            "Unexpected error running command: {0}".format("\n".join(e.cmd)))
    return result

def run_command(
    cmd, cmd_name, ignore_changed=False, final_target=False, config=False):
    if "--unlock" in cmd:
        # Bypass run command if --unlock is passed
        info("Unlocking working directory", color="blue")
        run_subprocess(cmd)
        return
    
    info("Running MTSv {0}".format(cmd_name))
    if ignore_changed:
        changed_targets = []
        warn("""
            Ignoring parameter changes.
            This is may cause down stream analysis to be wrong.
            Do not ignore parameter changes if binning parameters,
            kmer size, or sequence database have been modified.
            """)
    else:
        # try to capture snakemake utility options
        if set(PASSTHROUGH_ARGS).intersection(set(cmd)):
            sp.run(cmd)
            return
        info(
            "Checking if parameters have changed.",
            color="blue")
        changed_targets = get_targets_with_changed_params(cmd)
        if changed_targets:
            warn(
                "Parameters changed for targets: {0}\nRerunning".format(
                    "\n".join(changed_targets)))
        else:
            info("No parameters have changed.", color="blue")
    # add filter_candidate_taxa because it is not revaluated after 
    # checkpoint is hit.
    if final_target:
        changed_targets.append(final_target)
    cmd = cmd if not changed_targets else add_force_targets(
        cmd, changed_targets)
    if config:
        cmd = add_config(cmd, config)
    dryrun_flag = set(["--dryrun", "--dry-run", "-n"]).intersection(set(cmd))
    if dryrun_flag:
        # don't capture standard out for dryrun
        sp.run(cmd)
    else:
        p = run_subprocess(cmd)
        info("Finished running MTSv {0}".format(cmd_name))
    
def add_force_targets(cmd, changed_targets):
    """
    Add forced targets to command by appending to an existing
    --forcerun/-R argument or appending to the end of the command.
    cmd (list of strings): passed snakemake commands.
    changed_targets: list of targets that have parameter changes.
    """
    force_flag = list(set(["--forcerun", "-R"]).intersection(set(cmd)))
    if force_flag:
        idx = max([cmd.index(c) for c in force_flag]) # index of last occurance
        return cmd[: idx + 1] + changed_targets + cmd[idx + 1: ]
    else:
        # add forced targets to end if arg not already used.
        return cmd + ["-R"] + changed_targets

        
def add_config(cmd, config):
    """
    Adds config information to cmd either by appending to an existing
    --config/-C argument or appending to end of command.
    cmd (list of strings): passed snakemake commands.
    config (list of strings): additional values to add.
    """
    config_flag = list(set(["--config", "-C"]).intersection(set(cmd)))
    if config_flag:
        # append to end of config arg if passed
        idx = max([cmd.index(c) for c in config_flag]) # last occurance
        return cmd[: idx + 1] + config + cmd[idx + 1: ]

    else:
        # add config to end if arg not already used
        return cmd + ["--config"] + config 




def get_targets_with_changed_params(cmd):
    """
    Return a list of targets that have targets effected by parameter
    changes.
    """
    dryrun_flag = set(["--dryrun", "--dry-run", "-n"]).intersection(set(cmd))
    if dryrun_flag:
        # Dryrun will cause stdout to be inserted
        cmd = cmd[:]
        for flag in dryrun_flag:
            cmd.remove(flag)
        
    result = run_subprocess(cmd + ["--list-params-changes"])
    return [
        target for target in 
        result.stdout.decode('utf-8').strip().split("\n")
        if target]


def get_cmd(
    configfile_name, snakemake_args, until=False):
    cmd = [
            "snakemake", "--snakefile", SNAKEPATH]
    if until:
        cmd += ["--until", until]
    cmd += list(snakemake_args)
    cmd += ["--configfile", configfile_name]
    return cmd


def format_description(description):
    return "\n#".join(description.strip().split("  "))

def format_config_string(db_path):
    with open(SPECFILE_PATH, 'r') as handle:
        config = yaml.safe_load(handle)
    group_config = {
        "{}_description".format(k): format_description(v['description']) for
        k, v in config['groups'].items()}
    param_desc_config = {
        "{}_description".format(k): format_description(v['description'])
        for k, v in config['properties'].items()}
    param_default_config = {
        "{}_default".format(k): v.get('default', None)
        for k, v in config['properties'].items()}
    return CONFIG_STRING.format(
        database_path=db_path, datastore=data_path(),
        **group_config, **param_desc_config, **param_default_config)

@click.group(cls=OrderedGroup)
@click.version_option(
    version=VERSION, prog_name=PROG_NAME,
    message='%(prog)s: v%(version)s')
def cli():
    pass


@cli.command()
@click.option(
    '--configfile', '-c', type=click.File(mode='w'),
    default="mtsv.cfg",
    help='Specify path to write config file.', show_default=True)
@click.argument(
    'db_path', type=click.Path(exists=True, dir_okay=False, resolve_path=True))
def init(configfile, db_path):
    """
    Initialize new config file with default parameters.
    Requires path to sequence database config file.
    """
    configfile.write(format_config_string(db_path))
    click.secho(
        "MTSv config written to {0}".format(configfile.name),
        fg='green')


@cli.command()
@click.option(
    '--configfile', '-c', type=click.File(mode='w'),
    default="cluster.cfg",
    help='Specify path to write cluster config file.', show_default=True)
def cluster_init(configfile):
    """
    Initialize new cluster config file.
    Sections for each rule are provided with defaults but values
    should be modified depending on the size of the dataset.
    """
    configfile.write(CLUSTER_CONFIG)
    click.secho(
        "MTSv cluster config written to {0}".format(configfile.name),
        fg='green')



@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters."
)
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def pipeline(configfile, ignore_changes, snakemake_args):
    """
    Run MTSv pipeline from readprep to analyze.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE'
    """
    cmd = get_cmd(configfile.name, snakemake_args)
    run_command(cmd, "Pipeline", ignore_changed=ignore_changes)



@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
    )
def readprep(configfile, ignore_changes, snakemake_args):
    """
    QC reads, generate and deduplicate query kmers.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE --until readprep'
    """
    cmd = get_cmd(configfile.name, snakemake_args, until="readprep")
    run_command(cmd, "Readprep", ignore_changed=ignore_changes)


@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def binning(configfile, ignore_changes, snakemake_args):
    """
    Alignment-based metagenomic binning.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE --until collapse'
    """
    cmd = get_cmd(configfile.name, snakemake_args, until="collapse")
    run_command(cmd, "Binning", ignore_changed=ignore_changes)



@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters."
)
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def summary(configfile, ignore_changes, snakemake_args):
    """
    Summarize taxa hits.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE --until summary'
    """
    cmd = get_cmd(configfile.name, snakemake_args, until="summary")
    run_command(cmd, "Summary", ignore_changed=ignore_changes)



@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def analyze(configfile, ignore_changes, snakemake_args):
    """
    Perform statistical analysis on candidate taxa.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE'
    """
    cmd = get_cmd(configfile.name, snakemake_args)
    run_command(cmd, "Analyze", ignore_changed=ignore_changes)


@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--report', type=click.Path(exists=False, dir_okay=False, writable=True),
    default="results/report.html",
    help="Path to write report"
)
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def report(configfile, report, snakemake_args):
    """
    Generate report.\n

    Additional Snakemake parameters should be passed
    at the command line.
    
    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE --report results/report.html'
    """
    cmd = get_cmd(configfile.name,
        ["--report", report] + list(snakemake_args))
    run_command(
        cmd, "Report")
    info("")
    info("Report written to {}".format(report))


@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def extract(configfile, ignore_changes, snakemake_args):
    """
    Extract queries that hit a given taxid.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE -R extract
    --config run_extract=True CONFIG'
    """
    cmd = get_cmd(configfile.name, snakemake_args)
    run_command(
        cmd, "Extract",
        config=["run_extract=True"],
        ignore_changed=ignore_changes)

@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")    
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def extract_unaligned(configfile, ignore_changes, snakemake_args):
    """
    Extract queries that did not align to any taxa.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE -R unaligned_queries
    --config run_extract_unaligned=True CONFIG'
    """
    cmd = get_cmd(configfile.name, snakemake_args)
    run_command(cmd, "Extract Unaligned",
                config=["run_extract_unaligned=True"],
                ignore_changed=ignore_changes)

@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def wgfast(configfile, ignore_changes, snakemake_args):
    """
    SNP typing for strain-level resolution.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE -R wgfast
    --config run_wgfast=True CONFIG'
    """
    cmd = get_cmd(configfile.name, snakemake_args)
    run_command(cmd, "WGFAST", ignore_changed=ignore_changes,
                config=["run_wgfast=True"])


@cli.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option(
    '--configfile', '-c', type=click.File(mode='r'),
    required=True,
    help="Specify path to write config file.")
@click.option(
    '--ignore-changes', is_flag=True, default=False,
    help="Do not check for changes to config parameters.")
@click.argument(
    'snakemake_args', nargs=-1, type=click.UNPROCESSED,
)
def concoct(configfile, ignore_changes, snakemake_args):
    """
    Alignment-free binning of unaligned queries.\n
    Additional Snakemake parameters should be passed
    at the command line.

    Runs:\n
    'snakemake --snakefile SNAKEPATH
    --configfile CONFIGFILE -R fasta_bins
    --config run_concoct=True CONFIG'
    """
    cmd = get_cmd(configfile.name, snakemake_args)
    run_command(cmd, "Concoct",
                ignore_changed=ignore_changes,
                config=["run_concoct=True"])


@cli.command()
@click.argument(
    "datastore", type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.option(
    '--listing', '-ls', is_flag=True,
    help="List the datasets in the datastore and approximate size")
@click.option(
    '--remove', '-rm', type=click.STRING,
    help="Remove dataset from datastore (pass ID shown using --listing option)"
)
def datastore(datastore, listing, remove):
    """
    Manage datastore of expected value calculations.
    """
    if listing:
        
        table = get_datastore_info_table(datastore)

        click.secho(table.to_string(), fg="blue")
        # except:
        #     error("Cannot find datastore: {0}".format(datastore))
    elif remove is not None:
        warn(
            """
            The MTSv datastores keep previously collected expected values for 
            faster statistical analysis. Deleting datastores will require
            these to be recalculated. Datastores should be deleted if a new
            sequence database has been build and the expected values
            calculated from the previous database are no longer needed.
            """)

        while True:
            click.secho(
                'Are you sure you want to remove? [y/n] ', nl=False, fg="red")
            c = click.getchar()
            click.echo()
            if c.lower() == 'y':
                remove_datastore(datastore, remove)
                click.secho('Deleting: {0}'.format(remove), fg="red")
                break
            elif c.lower() == 'n':
                click.secho('Abort!', fg="red")
                break
            else:
                click.secho('Invalid input', fg="red")
        


if __name__ == '__main__':
    cli()
    
