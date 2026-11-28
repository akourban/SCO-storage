# Optimisation du stockage de donnée d'un code Monte-Carlo pour le calcul de la densité d'état à transition de spin

## Composition du groupe
### Etudiants , 
+ Santiago Vargas Triana, 
+ Anissa Kourban Houssen Fassy  
### Encadrant  
+ Thomas Dufaud

## Description

On s’intéresse à la simulation numérique des nano-matériaux à transition de spin (SCO), et plus particulièrement au calcul de la densité d’état de tels matériaux par une méthode de Monte Carlo.

Le projet consiste dans un premier temps à porter une solution de stockage de donnée développée pour le code Dyprosco dans un code MCES pour la simulation de matériaux 2D.  Puis on proposera une solution pour le stockage de tableaux de densité d’état dans le cas de matériaux 3D. Au second semestre on adaptera ces solutions à la version parallèle du code MCES. 

## Compilation

Depuis la racine du projet :

```bash
mkdir build
cmake -B build
cmake --build build
```

## Exécution

Depuis la racine du projet :
```bash
./build/mces
```

## Avancement au 28 novembre 2026

- Début du nettoyage et de la réorganisation du code source.
- Création du dossier `include` et du fichier d’en-tête `code_source.h`.
- Séparation de `main()` dans `src/main.c`.
- Mise en place de CMake pour compiler le projet.
- Passage des tableaux de densité et de l’histogramme en `double` (anciennement float).

### À poursuivre
- Tableaux `d` et `h` encore tridimensionnels ; les fonctions size2D et index2D ne sont pas encore utilisées pour les remplacer.
- Continuer le nettoyage et vérifier la compilation après chaque modification.
- Commencer le rapport (à rendre pour le 18/12)
- 
