Article MCES 2003 : 

Rappel : Ising-like model : 
on a des molécules contenant 3d⁴-3d⁷ ions de métal en symétrie octahédrique qui possèdent 2 états : 
(HS) : état high-spin de dégénérescence gHS
(LS) : état low-spin de dégénérescence gLS
(HS) et (LS) en fonction du champ de ligand, variations en fonction de la température, pression, champ magnétique...
HS -> 1, LS-> 0 (-1 si on est en physique mais on est en informatique), représentées par un spin σ
Hamiltonien qui contient les intéractions long-range(G)/short-range(J) écrit comme suit : 
H = -h Σ(σi)-J Σ(σiσj), 
avec h = -Δ/2 + (k[B]Tln(gHS/gLS))/2 + G<σ>
et Δ>0 ecart d'énergie du champ de ligand tq LS est l'état fondamental
fraction des moléculs en HS écrit comme suit : nHS = (1+<σ>)/2
Ising hamiltonien, calcul exact pour 1D

MCES : 
Energie totale du système : E(m,s) = -hm-Js 
avec m = Σ(σi), s =Σ(σiσj) 
MCES sert à calculer la densité d'état d(m,s): nombre de configurations pour un jeu de données de m et s
d(m,s) : dégénérescence de m et s
L'entropic sampling permet à n'importe quelle distribution P d'ếtre dérivée par un méthode de Monte-Carlo juste en introduisant la même distribution en tant que biais, en gros : 
PiW(i->j) = PjW(j->i)

Pour parcourir l'espace de tous les états, la méthode introduite juste avant doit être appliquée de la manière suivante : favoriser les densités d'état faibles et réduire les densités d'état fortes. Appliquer la probabilité des biais, cad inverse des densités d'état souhaitées : 
Pi∝1/(d(mi,si)), par conséquent : 
W(i->j)/w(j->i) = Pj/Pi =  (d(mj,sj))/(d(mi,si)) ou d(m,s) est à priori inconnue. 

Méthode interactive : résultats redimensionnés d'un système similaire à une échelle plus petite. Soit di(m,s) densité d'état après i itérations, si di(m,s) est utilisé en tant que biais on fait un sampling appelé 'Monte-Carlo stage' et renvoie un histogramme de la fréquence appelé H(m,s).
H(m,s)∝d(m,s)(1/di(m,s)). Après lacorrection avec le biais le resultat des densités d'état restreint est : di+1 (m,s)∝  di(m,s)Hi(m,s). 
Méthode itérative lancée qu'une fois pour remplir le tableau de H(m,s). Pour les paramètres suivants : température, interactions, écart d'énergie, dégénerescence, fonction Z = Σd(m,s) exp(-b(-hm-Js)) et toutes les propriétés thermodynamiques du système : (<σ>, capacité thermique, susceptibilité magnétique) peuvent être dérivées de manière analytique.
