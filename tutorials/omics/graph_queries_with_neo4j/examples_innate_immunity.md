# Innate Immunity Signaling — Example Queries

Exploring the TLR4 signaling pathway and its downstream transcription factors: NF-kB, JUN, and IRF3.

## Key proteins

| Protein | entryNameID   | Role                                        |
| ------- | ------------- | ------------------------------------------- |
| TLR4    | `TLR4_HUMAN`  | Toll-like receptor 4, pathogen recognition  |
| NF-kB1  | `NFKB1_HUMAN` | Transcription factor, inflammatory response |
| RelA    | `RELA_HUMAN`  | NF-kB subunit p65                           |
| JUN     | `JUN_HUMAN`   | AP-1 transcription factor component         |
| IRF3    | `IRF3_HUMAN`  | Interferon regulatory factor 3              |
| MYD88   | `MYD88_HUMAN` | TLR adaptor protein                         |
| TRAF6   | `TRAF6_HUMAN` | E3 ubiquitin ligase, NF-kB activation       |
| IRAK4   | `IRAK4_HUMAN` | IL-1 receptor-associated kinase 4           |

---

## Interaction Queries

### 1.1 What does TLR4 interact with?

```cypher
MATCH (p:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

or to view a graph, you can use `RETURN *` instead:

```cypher
MATCH (p:Protein {entryNameID: 'TLR4_HUMAN'})-[r:INTERACTS_WITH]-(partner:Protein)
RETURN *
```

### 1.2 What does NF-kB1 interact with?

```cypher
MATCH (p:Protein {entryNameID: 'NFKB1_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

### 1.3 Shortest path from TLR4 to NF-kB1

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'NFKB1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

or:

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'NFKB1_HUMAN'})
)
RETURN *
```

### 1.4 All shortest paths from TLR4 to NF-kB1

```cypher
MATCH path = allShortestPaths(
  (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'NFKB1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
ORDER BY hops
```

or

```cypher
MATCH path = allShortestPaths(
  (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[r:INTERACTS_WITH*..5]-(p2:Protein {entryNameID:
  'NFKB1_HUMAN'})
)
RETURN *
```

### 1.5 Shortest path from TLR4 to IRF3

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'IRF3_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 1.6 Shortest path from TLR4 to JUN

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..6]-(p2:Protein {entryNameID: 'JUN_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 1.7 Most connected proteins in the TLR4 subgraph

```cypher
MATCH (p:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..3]-(partner:Protein)
WITH partner, count(*) AS connections
RETURN partner.entryNameID, connections
ORDER BY connections DESC
LIMIT 10
```

### 1.8 Explore the TLR4 local subgraph

Best run in Neo4j Browser for graph visualization.

```cypher
MATCH path = (p:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..3]-(partner:Protein)
RETURN path
```

This will display a lot of nodes.

---

## Disease Queries

### 2.1 Diseases associated with TLR4

```cypher
MATCH (p:Protein {entryNameID: 'TLR4_HUMAN'})-[:ASSOCIATED_WITH]->(d:Disease)
RETURN d.name, d.id
```

=> No disease association

### 2.2 Diseases associated with NF-kB pathway proteins

```cypher
MATCH (p:Protein)-[:ASSOCIATED_WITH]->(d:Disease)
WHERE p.entryNameID IN ['TLR4_HUMAN', 'NFKB1_HUMAN', 'RELA_HUMAN', 'MYD88_HUMAN', 'IRAK4_HUMAN', 'TRAF6_HUMAN']
RETURN p.entryNameID, collect(d.name) AS diseases
```

=> Some diseases.

### 2.3 Cancer-associated proteins interacting with NF-kB1

```cypher
MATCH (nfkb:Protein {entryNameID: 'NFKB1_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)-[:ASSOCIATED_WITH]->(d:Disease)
WHERE d.name CONTAINS 'cancer'
   OR d.name CONTAINS 'carcinoma'
   OR d.name CONTAINS 'tumor'
RETURN DISTINCT partner.entryNameID, partner.proteinName, collect(d.name) AS diseases
```

### 2.4 Disease proteins that are hubs in the TLR4 network

```cypher
MATCH (tlr4:Protein {entryNameID: 'TLR4_HUMAN'})-[:INTERACTS_WITH*..3]-(p:Protein)-[:ASSOCIATED_WITH]->(d:Disease)
WITH p, collect(DISTINCT d.name) AS diseases, count(DISTINCT d) AS diseaseCount
MATCH (p)-[:INTERACTS_WITH]-(partner:Protein)
RETURN p.entryNameID, p.proteinName, count(partner) AS degree, diseaseCount, diseases
ORDER BY degree DESC
LIMIT 10
```

---

## GO Queries

### 3.1 Biological processes of TLR4

```cypher
MATCH (p:Protein {entryNameID: 'TLR4_HUMAN'})-[:INVOLVED_IN]->(bp:BiologicalProcess)
RETURN bp.name, bp.id
```

### 3.2 Shared biological processes between TLR4 and IRF3

```cypher
MATCH (p1:Protein {entryNameID: 'TLR4_HUMAN'})-[:INVOLVED_IN]->(bp:BiologicalProcess)<-[:INVOLVED_IN]-(p2:Protein {entryNameID: 'IRF3_HUMAN'})
RETURN bp.name, bp.id
```

### 3.3 Molecular functions of the NF-kB pathway proteins

```cypher
MATCH (p:Protein)-[:HAS_FUNCTION]->(mf:MolecularFunction)
WHERE p.entryNameID IN ['NFKB1_HUMAN', 'RELA_HUMAN', 'JUN_HUMAN', 'IRF3_HUMAN']
RETURN p.entryNameID, collect(mf.name) AS functions
```

### 3.4 Cellular components — where are TLR4 pathway proteins located?

```cypher
MATCH (p:Protein)-[:IS_LOCATED_IN]->(cc:CellularComponent)
WHERE p.entryNameID IN ['TLR4_HUMAN', 'NFKB1_HUMAN', 'RELA_HUMAN', 'MYD88_HUMAN', 'IRF3_HUMAN', 'JUN_HUMAN']
RETURN p.entryNameID, collect(cc.name) AS locations
```

or

```cypher
MATCH (p:Protein)-[r:IS_LOCATED_IN]->(cc:CellularComponent)
WHERE p.entryNameID IN ['TLR4_HUMAN', 'NFKB1_HUMAN', 'RELA_HUMAN', 'MYD88_HUMAN', 'IRF3_HUMAN', 'JUN_HUMAN']
RETURN *
```
