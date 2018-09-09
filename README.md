# monroe_bokma_2017

Code from Monroe, Melanie J., and Folmer Bokma. "Does density-dependent diversification mirror ecological competitive exclusion?." PloS one 12.10 (2017): e0184814.

Email from Bokma:

> hier [...] de matlab-code (testRho.m). Er hoort een ander scripts bij (Evo.m), en ik stuur ook wat data mee zodat je de code daar eerst eens op los kunt laten om te testen: Eerst haal je de soortnamen en lichaamsgrootte uit het bestand "dendroica_masses.csv"; de soortnamen stop je in een cell array, en de groottedata in een column vector. Dan importeer je de fylogenie: T=phytreeread('dendroica.nwk')
>
> En daarna start je de analyse: [ns,robs,p]=testRho(T,Dendroica_names,Dendroica_masses);
Het kan een tijdle duren want de code zoekt eerst de maximale correlatie voor de echte data, en daarna die voor 1000 simulaties op dezelfde boom. Als dat te lang duurt kun je in de code (testRho regel 16) de 1000 aanpassen. De output:
> ns: het aantal soorten
> robs: de maximale correlatie tussen de posities van de soorten in de fylogenie en hun (in dit geval) lichaamsgrootte
> p: de kans dat je een correlatie zo groot als robs krijgt als evolution willekeurig is met betrekking tot (in dit geval) lichaamsgrootte.
