
cd ../Data/
cat annot_annot.txt annot_size_leaf.txt annot_color_cluster.txt > ANNOT.txt
graphlan_annotate.py --annot ANNOT.txt tree.txt tree_bioreactor.xml
graphlan.py tree_bioreactor.xml tree_bioreactor.png --dpi 600 --size 10
