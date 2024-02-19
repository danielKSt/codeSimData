Vedlagt er en hdf5 fil som inneholder to arrays: 
  (1) Et array 'surf' som inneholder informasjon om virvler på overflaten
  (2) Et array 'sdiv' som inneholder informasjon om overflatedivergensen

(1) 'surf'
  Type: int8
  Størrelse: 12500x256x256
  Forklaring: Arrayverdier representerer en 256x256 overflate over 12500 tidssteg.
              Hver verdi er enten '0' eller '1': '1' representerer at det aktuelle
              gridpunktet er en virvel (eller del av en virvel), '0' representer
              at det ikke er en virvel der.

(2) 'sdiv'
  Type: single
  Størrelse: 12500x256x256
  Forklaring: Arrayverdier representerer en 256x256 overflate over 12500 tidssteg.
              Verdien i ett gridpunkt uttrykker overflatedivergensen, det vil si
              du/dx + dv/dy, i dette gridpunktet på det aktuelle tidspunktet.

Filen skal kunne lese i ditt foretrukne programmeringsspråk. Eksempelvis i Matlab
vil innlasting av data se noe slik ut:

    overflate_virvler = h5read("surface_vortices.h5", "/surf");
    overflate_divergens = h5read("surface_vortices.h5", "/sdiv");

Enn så lenge har vi bare bruk midlede verdier av overflate divergens i andre, og 
funnet en sammenheng mellom antall virvler (av en viss størrelse) som er tilstede
på overflaten og snitt overflate divergens per tidssteg. 

Fra dataen som er tatt med her vil dette være omtrent som å finne antall sammenkoblede
strukturer i arrayet 'surf' på et tidspunkt t, og plotte denne verdien mot verdien
av mean(s_div.^2,[2,3]), jf. Matlab-notasjon. 

Noe som skiller seg ut her er imidlertid at når vi detekterer virvler, filtrerer vi 
bort virvler som "lever" svært kort eller som er veldig små. I arrayet 'surf' er
alle virvler med, også dem som er veldig små og som lever svært kort.


Annet: 
- I dataen som er vedlagt nå er det ikke implementert noen tracking av virvler. 
  Hvert tidssteg gir altså bare informasjon om inneværende virvler, ikke hvordan
  disse er koblet til virvler på tidligere tidssteg. Vi har også slik informasjon
  tilgjengelig, og kan samle den og sende den om ønskelig.
- Noen virlver har en sentroide som ikke selv er i virvelen. Dette kan forekomme
  blant annet når to virvler kombineres til en større virvel, eller når virvelen
  har et "hull" i midten. Sistnevnte er en artefakt av hvordan vi setter "grensen" 
  for hva som skal defineres som en virvel (i literaturen er det flere konkurrerende 
  definisjoner). Rent fysisk vil virvelen koblet til overflaten ikke ha hull i midten.
  Skaper dette problemer for deres beregninger kan vi nok gjøre noe med denne dataen.


Ved spørsmål, ta kontakt: jorgen.r.aarnes@ntnu.no


