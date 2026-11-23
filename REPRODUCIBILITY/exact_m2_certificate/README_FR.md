# P030 — Certificat numérique BKM pour le bain orthogonal m=2

Pour `lambda=0.05`, `c=5`, `N=8500`, `gamma_s=1`, `gamma_f=3`,
`T=2000 ln(20)`, avec `C0=0.2 I`, `b=(0,0.1)` et `a0=0.6`,
les démonstrations et programmes de certification donnent

`3.95299e-41 < Delta_BKM,continuum < 3.98513e-41`.

Le résultat est **assisté par ordinateur** et limité à ce point. Les calculs
reposent sur IEEE 754 binary64 (arrondi au plus proche, sans contraction FMA
implicite ni réassociation), sur les bornes explicites et les bibliothèques
de calcul exact/par intervalles identifiées dans les méthodes.

Lire `README_REPRODUCE.md`, `CERTIFICATE_IDEAL_ORDER2_2026-11-22.md`,
`PROOF_DD_AND_NUMERICAL_CERTIFICATE_2026-11-22.md` et
`continuum_proof/CONTINUUM_BOUND_AND_TRANSFER.md`.
Les deux fichiers d'entrée sont inclus et contrôlés par SHA-256.
Les calculs anciens obtenus avec les coefficients SVD arrondis ne sont pas
couverts par le certificat du polynôme exactement orthogonal.
