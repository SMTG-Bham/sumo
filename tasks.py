"""

Deployment file to facilitate releases of sumo.
Note that this file is meant to be run from the root directory of the repo.
"""

from invoke import task

import os
import json
import requests
import re

__author__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Oct 20 2017"


@task
def publish(ctx):
    ctx.run("rm dist/*.*", warn=True)
    ctx.run("python3 setup.py sdist bdist_wheel")
    ctx.run("twine upload dist/*")


@task
def release(ctx):
    with open("CHANGELOG.rst") as f:
        contents = f.read()
    new_ver = re.findall('\n(v.*)', contents)[0]
    toks = re.finditer("------\n*(.*?)^\`\[v\d\.\d\.\d", contents,
                       re.MULTILINE | re.DOTALL)
    desc = list(toks)[1].groups()[0].strip()
    payload = {
        "tag_name": new_ver,
        "target_commitish": "master",
        "name": new_ver,
        "body": desc,
        "draft": False,
        "prerelease": False
    }
    response = requests.post(
        "https://api.github.com/repos/smtg-ucl/sumo/releases",
        data=json.dumps(payload),
        headers={"Authorization": "token " + os.environ["GITHUB_TOKEN"]})
    print(response.text)
