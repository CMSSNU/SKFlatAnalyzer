import os,sys,re
if __name__=="__main__":
    ForSNU=os.environ["DATA_DIR"]+"/*/Sample/ForSNU"
    samples=os.popen('find {} -maxdepth 1 -mindepth 1 -name "*.txt"'.format(ForSNU)).read().split()
    
    for sample in samples:
        samplename=os.path.basename(sample).replace(".txt","")
        era=re.search(r"/{tag}/([^/]*)/".format(tag=os.environ["SKFlatV"]),sample).group(1)
        last=int(os.popen("cat {} | sort -V |tail -n1|xargs -i basename {{}}|sed 's/^.*_\([0-9]*\)\.root$/\\1/'".format(sample)).read())
        num=int(os.popen("cat {} | wc -l".format(sample)).read())
        if last!=num:
            print era,samplename,":",num,"/",last
