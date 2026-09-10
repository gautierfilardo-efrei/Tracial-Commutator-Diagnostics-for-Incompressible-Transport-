# Créer le dépôt demandé

Compte demandé : `gautierfilardo-efrei`.
Nom proposé : `tracial-transport-diagnostics`.
Visibilité initiale du script : **privée**.

Cette archive ne constitue pas une création sur GitHub. Aucune connexion au
compte n'était active au moment de sa préparation.

Avec une connexion GitHub disponible dans ChatGPT, la création et l'envoi des
fichiers peuvent être effectués directement à partir de ce projet.

Sinon, sur votre ordinateur, installez Git et GitHub CLI, authentifiez-vous avec
`gh auth login --hostname github.com`, puis vérifiez le compte avec
`gh api user --jq .login`. Ne transmettez aucun mot de passe ou jeton dans une
conversation. Configurez votre identité Git avec votre nom et l'adresse que
vous souhaitez associer à vos commits, si ce n'est pas déjà fait.

Après extraction de l'archive, exécutez depuis sa racine :

```sh
bash create_repository.sh
```

Le script vérifie le compte actif et la configuration Git, refuse d'écraser un
dépôt existant, initialise le projet local puis crée le dépôt privé et envoie
le premier commit. Il n'a pas été exécuté contre GitHub pendant la préparation.
S'il échoue après la création distante, inspectez l'état du dépôt avant toute
nouvelle tentative; aucun envoi forcé n'est prévu.

La commande de création suit la
[documentation officielle de GitHub CLI](https://cli.github.com/manual/gh_repo_create).

Après création, vérifiez l'onglet Actions. Associez ensuite le manuscrit au
commit réellement utilisé. L'ouverture publique, le choix d'une licence et
l'attribution éventuelle d'un DOI doivent correspondre à des décisions et des
enregistrements effectifs.
