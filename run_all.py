from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parent/"src"))
from experiments import run

if __name__=="__main__":
    p=argparse.ArgumentParser(description="Thin-Si self-consistent Schrodinger-Poisson experiments")
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parent/"figures")
    p.add_argument("--spacing",type=float,default=.05,help="Target silicon mesh spacing in nm")
    a=p.parse_args();run(a.output,a.spacing)
