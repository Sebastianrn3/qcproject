Šiame aplankale yra dar env skriptų B (2) ir C (1) veikimui+pradiniai .xyz, nėra tik paties mopac.
Naudojimas:

Pirmas scriptas/funkcijos: rotoriai:
1. Globaliuose konstantuose paruošiamas (rankiniu būdu) listas grupių (rotorių) su kiekvieno atomų numeriais, to paties ilgio inkarų listas. Paverčiama į 0-based
2. prepare_rigid_rotors: Kiekvienai grupei formuojamas objektas su ją sudarančiais atomų numeriais, inkarinio atomo id ir jo koordinates, kiekvieno atomo koordinatės inkaro atžvilgiu
1SCF:
3. build_geom_from_rotors: Surenkama pilna geometrija (Rodrigo rotacijos f-lė)
4. mopac atlieka 1SCF
5. pack_gradients_multi_rotor: Grupių atomų gradientai perskaičiuojami į 3 dE/dw kiekvienam rotoriui vietoje 3N

-----------
Antras skriptas: interpoliacija su optimizavimu slepiant fiksuotus atomus:
1. Nuskaitomi starto ir galo .xyz, tiesiškai interpoliuojama į n intermediatų
2. Patikrinamas fiksuotų atomų atitikimas starte ir gale, formuomajas bool listas kurie atomai optimizuojami

-create_dataset_of_images_from_xyz_files_via_1scf(): skaičiuojami 1SCF kiekvienai iš n+2 struktrūrų ir įrašo į .npz
-create_dataset_of_images_from_xyz_files_via_optimizer(): tas pats plius relaksuoja į lokalų mininumą

Optimizatoriui paduodami tik laisvų atomų parametrai, po iteracijos surenkama pilna geometrija, siunčiama mopac, nukerpa fiksuotų atomų parametrus naujai iteracijai
-----------
Antras scriptas optimizavo R, 3 intemediate images, P per 40 minčių. 
Pirmo skripto tik patikrai (ne fizikinei prasmei) irgi perleidau per adaptuotą optimizatorių ir su 9 ridinėmis grupėmis (2 his, 2 asn, 4 acetatų, H2O: inkaruojami pirmi atomai) - optimizavimas užtruk


0) Paleidimui reikia tik nustatyti  interpol.py faile MOPAC_EXE_PATH - kelią iki mopac.exe. Su naujausia mopac versija turi pakakti, kad zipe esantis aplankas būtų toje pačioje direktorijoje, kur ir mopac aplankas. 
1) interpol.py gale yra trys funkcijos, prieš naudojimą iškomentuojama reikiama. Pirma sukuria N_IMAGES image'ų ir visoms paskaičiuoja E, grad ir su kita informacija įrašo į images_raw.npz aplanke runs, kurią sukurus galima terminale pasižiūrėti su trečia funkcija.
Antra funkcija kiekvieną image relaksuoja optimizatoriuje bei konsolėje rodo, arčiau kurios pusės nuvažiavo geometrija (ir kiek relaksuota geometrija atitinka reaktantą/produktą). Antros funkcijos nerekomenduoju leisti (ilgas skaičiavimas + labiau skirtas tikrinimui/lyginimui), bet jei vis reiktų paskaityti ir jos npz rezultatą - jis bus įraše images_raw_opt.npz po pirmo paleidimo.
2) Norint išmėginti rigid_bodies_v2.py pakanka jo tūrinį įterpti į interpol.py galą (pastarojo paskutinės trys funkcijos turi likti užkomentuotos arba pašalintos).
-------------