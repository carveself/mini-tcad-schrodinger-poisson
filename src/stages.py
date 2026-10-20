"""Run a single reproduction checkpoint."""
from pathlib import Path
import argparse,json
from experiments import baseline,parameters,verification,feedback

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("stage",choices=["baseline","geometry","verify","feedback"])
    p.add_argument("--output",type=Path,default=Path("reproduction/checkpoint"))
    p.add_argument("--spacing",type=float,default=.05)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    f={"baseline":baseline,"geometry":parameters,"verify":verification,"feedback":feedback}[a.stage]
    result=f(a.output,a.spacing)
    if a.stage in ["verify","feedback"]:
        (a.output/(a.stage+".json")).write_text(json.dumps(result,indent=2)+"\n")
    print(f"Completed {a.stage}: {a.output.resolve()}")
