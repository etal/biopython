#!/usr/bin/python
import os

print ''
print ''
print '==============================='
print 'Testing TreeSum module, which makes use of revised Bio.Nexus.Trees module'
print '==============================='


# Checking tree functions
# Using update bug fix function by Brad Chapman
from Bio.Nexus.Trees import Tree
from Bio.Geography.TreeSum import TreeSum
#from TreeSum import TreeSum

tree_str = '(((((((((((((((((Sambucus:43.136024,Viburnum:43.136040)Adoxaceae:53.892513,(Acanthopanax:34.719704,Aralia:34.719727,Dendropanax:34.719727,Evodiopanax:34.719727,Kalopanax:34.719727,Schefflera:34.719727)Araliaceae:62.308830):7.045975,Ilex:104.074516):3.056864,((((((Catalpa:22.623766,Paulownia:22.623785)Bignoniaceae:22.623766,(Clerodendrum:19.864199,Premna:19.864218)Verbenaceae:25.383331):22.378326,(Chionanthus:29.443968,Forestiera:29.443979,Fraxinus:29.443979,Ligustrum:29.443979,Osmanthus:29.443979,Syringa:29.443979)Oleaceae:38.181892):19.113832,(Adina:38.252457,Cephalanthus:38.252472,Emmenopterys:38.252472,Pinckneya:38.252472,Randia:38.252472)Rubiaceae:48.487236):2.360018,Ehretia:89.099709):13.495450,Eucommia:102.595161):4.536214):0.905059,((((Clethra:78.134140,((Cliftonia:38.402752,Cyrilla:38.402775)Cyrillaceae:38.402752,(Arbutus:38.402752,Elliottia:38.402775,Enkianthus:38.402775,Kalmia:38.402775,Lyonia:38.402775,Oxydendrum:38.402775,Rhododendron:38.402775,Vaccinium:38.402775)Ericaceae:38.402752):1.328631):12.980787,(((Halesia:30.391993,Pterostyrax:30.392012,Styrax:30.392012)Styracaceae:51.775261,Symplocos:82.167252):0.000000,(Camellia:41.083626,Franklinia:41.083649,Gordonia:41.083649,Stewartia:41.083649,Ternstroemia:41.083649)Theaceae:41.083626):8.947675):0.000149,Diospyros:91.115099):2.023849,((Ardisia:18.344650,Myrsine:18.344666)Myrsinaceae:74.794174,Bumelia:93.138824):0.000101):14.897509):1.462594,((Alangium:48.167362,Aucuba:48.167370,Cornus:48.167370,Macrocarpium:48.167370,Torricellia:48.167370)Cornaceae:53.025345,(Hydrangea:97.032310,(Davidia:48.516151,Nyssa:48.516167)Nyssaceae:48.516151):4.160399):8.306321):7.064716,Schoepfia:116.563736):0.000000,((((Altingia:50.813206,Liquidambar:50.813213)Altingiaceae:50.813206,(Disanthus:50.813206,Distylium:50.813213,Fortuneria:50.813213,Hamamelis:50.813213,Loropetalum:50.813213,Sinowilsonia:50.813213)Hamamelidaceae:50.813206):0.000131,(Cercidiphyllum:87.828712,Daphniphyllum:87.828712):13.797829):13.247040,(((((((Choerospondias:21.440735,Cotinus:21.440742,Pistacia:21.440742,Rhus:21.440742,Toxicodendron:21.440742)Anacardiaceae:37.304596,(Acer:29.372665,Aesculus:29.372681,Dipteronia:29.372681,Koelreuteria:29.372681,Sapindus:29.372681)Sapindaceae:29.372665):0.000114,((Cedrela:49.350353,(Ailanthus:24.675177,Leitneria:24.675188,Picrasma:24.675188)Simaroubaceae:24.675177):4.016092,(Evodia:26.683222,Phellodendron:26.683233,Ptelea:26.683233,Zanthoxylum:26.683233)Rutaceae:26.683222):5.379002):29.842871,(Firmiana:32.917126,Tilia:32.917149)Malvaceae:55.671188):12.661992,(Lagerstroemia:84.110847,Szyzygium:84.110847):17.139463):2.612011,((((((Alnus:16.609535,Betula:16.609543,Carpinus:16.609543,Corylus:16.609543,Ostrya:16.609543)Betulaceae:37.306709,((Carya:25.504854,Cyclocarya:25.504866,Juglans:25.504866,Engelhardtia:25.504866,Platycarya:25.504866,Pterocarya:25.504866)Juglandaceae:25.504854,Myrica:51.009708):2.906531):9.893459,(Castanea:31.904850,Castanopsis:31.904873,Cyclobalanopsis:31.904873,Fagus:31.904873,Lithocarpus:31.904873,Quercus:31.904873)Fagaceae:31.904850):21.681023,(((((Celtis:20.739927,Pteroceltis:20.739939)Cannabaceae:20.739927,((Broussounetia:12.614990,Cudrania:12.615005,Maclura:12.615005,Morus:12.615005)Moraceae:12.614990,Oreocnide:25.229980):16.249876):10.909924,(Aphananthe:26.194889,Hemiptelea:26.194897,Planera:26.194897,Ulmus:26.194897,Zelkova:26.194897)Ulmaceae:26.194889):11.649286,(Hovenia:32.019470,Rhamnus:32.019493,Ziziphus:32.019493)Rhamnaceae:32.019596):8.938065,(((Amelanchier:36.488564,(Crataegus:36.488586,Mespilus:36.488586):0.000000):0.000000,Chaenomeles:36.488586,Eriobotrya:36.488586,Malus:36.488586,Photinia:36.488586,Pyrus:36.488586,Sorbus:36.488586):0.000000,Prunus:36.488586)Rosaceae:36.488564):12.513593):4.616908,(Albizia:31.901920,Cercis:31.901943,Cladrastis:31.901943,Dalbergia:31.901943,Erythrina:31.901943,Gleditsia:31.901943,Gymnocladus:31.901943,Laburnum:31.901943,Maackia:31.901943,Ormosia:31.901943,Robinia:31.901943,Sophora:31.901943)Fabaceae:58.205711):4.139401,((Euonymus:90.433327,Sloanea:90.433327):0.000101,((Mallotus:28.689901,Sapium:28.689920)Euphorbiaceae:50.330055,(Idesia:29.019764,Poliothyrsis:29.019779,Populus:29.019779,Salix:29.019779,Xylosma:29.019779)Salicaceae:50.000195):11.413469):3.813607):9.615288):0.000000,(Staphylea:21.372393,Tapiscia:21.372404,Turpinia:21.372404)Staphyleaceae:82.489929):11.011259):1.690163):7.829397,Buxus:124.393143):0.000000,Tetracentron:124.393143):2.763555,Meliosma:127.156693):1.664427,Platanus:128.821121):2.029122,Euptelea:130.850250):11.447736,((Asimina:95.972672,(Liriodendron:47.125092,Magnolia:47.125114,Manglieita:47.125114,Michelia:47.125114)Magnoliaceae:48.847580):46.325292,(Actinodaphne:49.903526,Cinnamomum:49.903542,Lindera:49.903542,Litsea:49.903542,Machilus:49.903542,Neolitsea:49.903542,Nothaphoebe:49.903542,Persea:49.903542,Phoebe:49.903542,Sassafras:49.903542,Umbellularia:49.903542)Lauraceae:92.394188):0.000257):1.840266,(Yucca:110.138222,((Sabal:100.000000,(Serenoa:95.000000,Trachycarpus:95.000000)ST:5.000000)Arecaceae:10.000000,(Arundinaria:20.476601,Phyllostachys:20.476624,Semiarundinaria:20.476624)Poaceae:89.661629):0.000000):34):30.861772,Illicium:175.000000)aus2ast:175.000000,(((((Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)TT1:25.000000)Taxaceae:90.000000,((((((((Calocedrus:85.000000,Platycladus:85.000000)CP:5.000000,(Cupressus:85.000000,Juniperus:85.000000)CJ:5.000000)CJCP:5.000000,Chamaecyparis:95.000000)CCJCP:5.000000,(Thuja:7.870000,Thujopsis:7.870000)TT2:92.13)CJCPTT:30.000000,((Cryptomeria:120.000000,Taxodium:120.000000)CT:5.000000,Glyptostrobus:125.000000)CTG:5.000000)CupCallTax:5.830000,((Metasequoia:125.000000,Sequoia:125.000000)MS:5.000000,Sequoiadendron:130.000000)Sequoioid:5.830000)STCC:49.060001,Taiwania:184.889999)Taw+others:15.110000,Cunninghamia:200.000000)nonSci:15.000000)Tax+nonSci:10.000000,Sciadopitys:225.000000):25.000000,(((Abies:106.000000,Keteleeria:106.000000)AK:54.000000,(Pseudolarix:156.000000,Tsuga:156.000000)NTP:4.000000)NTPAK:24.000000,((Larix:87.000000,Pseudotsuga:87.000000)LP:81.000000,(Picea:155.000000,Pinus:155.000000)PPC:13.000000)Pinoideae:16.000000)Pinaceae:66.000000)Coniferales:25.000000,Ginkgo:275.000000)gymnosperm:75.000000)seedplant:50.000000;'


tree_obj = Tree(tree_str)

print tree_obj
bigtree_ts = TreeSum(tree_obj)






# This works

trstr0 = "(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;"

to0 = Tree(trstr0)
print ''
print to0



# Gymnosperms tree with node labels; doesn't work
trstr1a = '(((((Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)TT1:25.000000)Taxaceae:90.000000,((((((((Calocedrus:85.000000,Platycladus:85.000000)CP:5.000000,(Cupressus:85.000000,Juniperus:85.000000)CJ:5.000000)CJCP:5.000000,Chamaecyparis:95.000000)CCJCP:5.000000,(Thuja:7.870000,Thujopsis:7.870000)TT2:92.13)CJCPTT:30.000000,((Cryptomeria:120.000000,Taxodium:120.000000)CT:5.000000,Glyptostrobus:125.000000)CTG:5.000000)CupCallTax:5.830000,((Metasequoia:125.000000,Sequoia:125.000000)MS:5.000000,Sequoiadendron:130.000000)Sequoioid:5.830000)STCC:49.060001,Taiwania:184.889999)Taw+others:15.110000,Cunninghamia:200.000000)nonSci:15.000000)Tax+nonSci:10.000000,Sciadopitys:225.000000):25.000000,(((Abies:106.000000,Keteleeria:106.000000)AK:54.000000,(Pseudolarix:156.000000,Tsuga:156.000000)NTP:4.000000)NTPAK:24.000000,((Larix:87.000000,Pseudotsuga:87.000000)LP:81.000000,(Picea:155.000000,Pinus:155.000000)PPC:13.000000)Pinoideae:16.000000)Pinaceae:66.000000)Coniferales:25.000000,Ginkgo:275.000000)gymnosperm:75.000000;'

to1a = Tree(trstr1a)
print ''
print to1a


# Just Taxaceae; doesn't work
trstr1b = '(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)TT1:25.000000)Taxaceae:90.000000;'
to1b = Tree(trstr1b)
print ''
print to1b

# Just Taxaceae; this works; node labels deleted
trstr1c = '(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)25.000000)90.000000;'
to1c = Tree(trstr1c)
print ''
print to1c



# This doesn't work (from bug report)
trstr2 = "(((t9:0.385832, (t8:0.445135,t4:0.41401)C:0.024032)B:0.041436, t6:0.392496)A:0.0291131, t2:0.497673, ((t0:0.301171, t7:0.482152)E:0.0268148, ((t5:0.0984167,t3:0.488578)G:0.0349662, t1:0.130208)F:0.0318288)D:0.0273876);"
to2 = Tree(trstr2) 
print ''
print to2








print ''
print ''
print '==============================='
print 'Testing TreeSum object, treefunctions.py, and revised Bio.Nexus.Trees module'
print '==============================='

ts = TreeSum(to2)


print ''
ts.test_Tree()


print ''
print 'Extracting subtree...'

# Node to make new root
newroot_id = ts.treeobj.chain[9].id

#newtree = copy.deepcopy(ts.treeobj)
#ts_newtree = TreeSum(newtree)
#ts_newtree.test_Tree()


list_leafids_to_keep = ts.get_successor_leafids(newroot_id)
print ''
print 'list_leafids_to_keep:', list_leafids_to_keep
newtree = ts.subset_Tree(list_leafids_to_keep)
ts_newtree = TreeSum(newtree)
ts_newtree.parent_tree = ts.treeobj
ts_newtree.test_Tree()


# This could be slow
#bigtree_ts.test_Tree()





