import os
from os.path import exists, join
from shutil import copy
import subprocess


def get_system_info():
    env = os.environ

    if exists("/home/jovyan/team_classify"):
        team_classify_dir = "/home/jovyan/team_classify"
        classifier = env["JUPYTERHUB_USER"]
        env["system"] = "tacc-cloud"
    elif exists("/work/05350/ecooper/wrangler/team_classify"):
        team_classify_dir = "/work/05350/ecooper/wrangler/team_classify"
        classifier = env["USER"]
        env["system"] = "tacc-work"
    elif exists("/home/idies/workspace/hetdex_workspace/team_classify"):
        team_classify_dir = "/home/idies/workspace/hetdex_workspace/team_classify"
        classifier = ""
        env["system"] = "sciserver"
    else:
        print("Can't determine your system. Manually enter classify directory")
        return None

    env["team_classify_dir"] = team_classify_dir
    env["classifier"] = classifier
    env["training_dir"] = join(team_classify_dir, "training")
    env["elixer_dir"] = join(team_classify_dir, "all_pngs")

    return env


def activate_file(filename, save_to_shared=True):
    """
    Script to copy over the classification file and
    either save to current working directory or to the shared
    space on Erin's /work path

    Parameters
    ----------
    filename
        filename str
    save_to_shared
        boolean flag to save to shared space. Default
        is True. If False it will save in the current
        working directory 
    """
    env = get_system_info()

    if (env["system"] == "tacc-cloud") and save_to_shared:
        your_classify_dir = join(
            env["team_classify_dir"], "classified", env["JUPYTERHUB_USER"]
        )
        if not exists(your_classify_dir):
            os.mkdir(your_classify_dir)
            subprocess.call(["chmod", "g+rw", your_classify_dir])
    elif (env["system"] == "tacc-work") and save_to_shared:
        your_classify_dir = join(env["team_classify_dir"], "classified", env["USER"])
    else:
        your_classify_dir = os.getcwd()

    env["your_classify_dir"] = your_classify_dir

    filepath = join(join(env["team_classify_dir"], "dets", filename))
    yourfile = join(your_classify_dir, filename)

    if not exists(yourfile):
        copy(filepath, yourfile)
        print("Copied over file to " + yourfile)
    else:
        print("File exists already at " + yourfile)

    return yourfile


def delete_file(yourfile):

    env = get_system_info()

    if exists(yourfile):
        print("Deleting " + yourfile)
        subprocess.call(["rm", yourfile])
    else:
        try:
            filename = yourfile
            yourfile = join(env["your_classify_dir"], filename)
            print("Deleting " + yourfile)
            subprocess.call(["rm", yourfile])
        except KeyError:
            print("Could not delete " + yourfile)
