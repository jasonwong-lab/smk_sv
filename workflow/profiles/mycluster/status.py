#!/usr/bin/env python3
"""
To get the status of a job in the CPOS PBS cluster.
"""

import sys
import subprocess

jobid = sys.argv[1]

try:
    res = subprocess.run(
        f"qstat -f -x {jobid}",
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
    )

    output = res.stdout.decode()
    lines = output.split('\n')
    data = {}
    for line in lines:
        parts = line.split(' = ')
        if len(parts) == 2:
            key, value = parts
            data[key.strip()] = value.strip()
    job_state = data.get("job_state", "")

    if job_state == "F":
        exit_status = data.get("Exit_status", "")
        if exit_status == "0":
            print("success")
        else:
            print("failed")
    else:
        print("running")

except (subprocess.CalledProcessError, IndexError, KeyboardInterrupt) as e:
    print("failed")
