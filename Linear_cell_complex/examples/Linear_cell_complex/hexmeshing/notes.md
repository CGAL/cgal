Notre implémentation:

1. Générer la grille
2. Identifier les cellules a rafiner selon un prédicat et un modèle chargé 
3. Boucle, pour chaque dim i, j, k 
* 3.1 Itérer sur chaque plan paire de la grille (xy si i, ....)
  * On itère sur chaque vertex du plan et on stocke pour chaque face, le nombre de noeuds marqués
  * On stocke les faces ayant été marqués au moins une fois 
  * On stocke également tous les vertices marqués 

 TODO Etape de propagation

* 3.2 On fait une deuxième passe mais cette fois ci sur les faces
  * On récupère toutes les faces incidentes hormis celles déjà sur la grille (càd qu'on va les traiter plus tard),
  Ces faces seront à traiter pour la substitution
  * On stocke également les deux volumes adjacents ou va se faire la substitution 
  

Puis tout ça en complexité lineaire 
* 3.3 Insertion des vertices 
* 3.4 Substitution des faces
* 3.5 Substitution des volumes 
    * 3.5.2 Gérer le cas spécifique du 3 template 


TODO: Dans les itérations suivantes varient légèrement
stocker le template utilisé dans un attribut de volume 
3.2 