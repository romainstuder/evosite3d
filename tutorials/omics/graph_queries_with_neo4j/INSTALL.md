# Installing Neo4j

## Option A — Neo4j Desktop (recommended for beginners)

1. Download Neo4j Desktop from [https://neo4j.com/download](https://neo4j.com/download)
2. Install and launch the application
3. Create a new project and add a local DBMS
4. Choose Neo4j version **5.x**, set a password, and start the database
5. Open **Neo4j Browser** from the interface to run Cypher queries

## Option B — Neo4j Community Server

```bash
# macOS with Homebrew
brew install neo4j
brew services start neo4j

# Linux (Ubuntu/Debian)
wget -O - https://debian.neo4j.com/neotechnology.gpg.key | sudo apt-key add -
echo 'deb https://debian.neo4j.com stable latest' | sudo tee /etc/apt/sources.list.d/neo4j.list
sudo apt update && sudo apt install neo4j
sudo systemctl start neo4j
```

Access Neo4j Browser at `http://localhost:7474` — default credentials are `neo4j / neo4j`.

## Option C — Docker

```bash
docker run \
  --name neo4j \
  -p 7474:7474 -p 7687:7687 \
  -e NEO4J_AUTH=neo4j/password \
  neo4j:latest
```

```bash
more /opt/homebrew/Cellar/neo4j/2026.03.1/libexec/conf/neo4j.conf
subl /opt/homebrew/Cellar/neo4j/2026.03.1/libexec/conf/neo4j.conf
```

```
=> Add/uncomment: initial.dbms.default_database=neo4j
=> Add/uncomment: dbms.security.allow_csv_import_from_file_urls=true
=> Add/uncomment: dbms.security.auth_enabled=false
=> Remove/Comment out: server.directories.import=import
```
