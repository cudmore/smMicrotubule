read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    # do dangerous stuff
    python deleteMaskAndSkel.py batch=../analysis/wt-female.txt
    python deleteMaskAndSkel.py batch=../analysis/wt-male.txt
    python deleteMaskAndSkel.py batch=../analysis/ko-female.txt
    python deleteMaskAndSkel.py batch=../analysis/ko-male.txt
    

fi

