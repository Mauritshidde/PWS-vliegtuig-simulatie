# vliegtuig Simulatie

## inhouds opgave
- [over het project](##over-het-project)
- [gebruikte library's](##gebruikte-library's)
- [build and installation](##build-and-installation)
- [uitleg over de controls](## controls)

## over het project
Voor school moesten we een PWS (profielwerkstuk) maken, dit PWS moest te maken hebben met 1 of meer van de vakken die we volgen. 
Het maken van de simulatie en het verslag wat we erover moesten maken heeft ons meer dan 80 uur per persoon gekost.
Tijdens de pws-dagen wisten we al vrij snel dat we de vakken Informatica/Coderclass en Natuurkunde wilde combineren in de vorm van een simulatie. 
Wij, Maurits en Kaio, hebben allebei interesse in het samenvoegen van deze vakken omdat we het leuk vinden om de Natuurkunde op een praktische manier te gebruiken. 
Het maken van een simulatie is een goede manier om een sterke intuïtie op te bouwen voor de natuurkundige concepten, en ook een goede uitdaging in op het gebied van de Informatica. 
We vonden het allebei een goed idee om de simulatie over een vliegtuig te laten gaan, omdat we bij natuurkunde nog niet veel over dit onderwerp geleerd hadden, ondanks dat we er wel interesse in hadden. 
Ook biedt het onderwerp de mogelijkheid om de simulatie te verkleinen of complexer te maken, als het blijkt dat onze originele doelen niet helemaal voldoen. 
Bij het schrijven van code is het vaak moeilijk om in te schatten hoelang iets zal duren, omdat onverwachte errors of moeilijkheden kunnen opduiken.  

Uiteindelijk hebben we besloten op de volgende opzet:   

We maken een simulatie van een vliegtuig, in de programmeertaal C++, waarbij het vliegtuig in stabiele vlucht is.
Tijdens het uitvoeren van de simulatie kun je een aantal variabelen in real-time aanpassen, zoals de kracht van de motor. 
Voor de luchtweerstand en de liftkracht van het vliegtuig berekenen we nodige constanten van tevoren in een apart programma. 

## gebruikte library's
We hebben in dit project een aantal library's gemaakt, voor het visualiseren en plotten.

#### raylib
We hebben gerbuik gemaakt van raylib voor het visualiseren van de simulatie. op de github pagina van raylib staat hoe je het kan instaleren. Hieronder staat ook een van de mogelijke manieren voor het instaleren van raylib. Het project maakt ook gerbuik van raygui en raymath dit zijn sub librarys van raylib en zijn al in het project gestopt.
##### instalatie voor de static versie:
Je kan raylib instaleren met de volgende commands of mogelijk via de package manager van je linux distro.
```
git clone https://github.com/raysan5/raylib.git raylib
cd raylib/src/
make PLATFORM=PLATFORM_DESKTOP
```
##### de link naar de website van raylib is:
https://www.raylib.com/
##### de link naar de github page van raylib is:
https://github.com/raysan5/raylib

#### matplotlib-cpp
Het project maakt gerbuik van matplotlib-cpp, voor het plotten van grafieken over het vleigtuig.

##### instalatie
```
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
./vcpkg integrate install
vcpkg install matplotlib-cpp
```

#### nlohmann json
Voor het opslaan van gegevens over een vliegtuig en het verkrijgen van de lift wordt er gebruik gemaakt van json. Hiervoor is de libray nlohmann json nodig. Je kan de library instaleren via de package manager van de linux distrubutie die je gebruikt.


##### build and installation

Voor het instaleren van het project moet je het project downloaden en unzippen.   
Voor het builden moet je in de folder van het project zitten en dan ```make``` uitvoeren in de terminal.
Voor het uitvoeren van het programma moet je dan ```./a.out``` uitvoeren in de terminal.

## running


