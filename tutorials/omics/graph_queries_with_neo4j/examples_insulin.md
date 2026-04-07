# Insulin Signaling Pathway — Example Queries

The insulin signaling pathway is a textbook example well suited for graph exploration.

## Key proteins

| Protein          | entryNameID    | Role                           |
| ---------------- | -------------- | ------------------------------ |
| Insulin receptor | `INSR_HUMAN`   | Receptor tyrosine kinase       |
| IRS1             | `IRS1_HUMAN`   | Docking/adaptor protein        |
| PI3K regulatory  | `PIK3R1_HUMAN` | Lipid kinase complex           |
| PI3K catalytic   | `PIK3CA_HUMAN` | Lipid kinase complex           |
| AKT1             | `AKT1_HUMAN`   | Serine/threonine kinase        |
| mTOR             | `MTOR_HUMAN`   | Growth & metabolism regulator  |
| GLUT4            | `SLC2A4_HUMAN` | Glucose transporter            |
| GRB2             | `GRB2_HUMAN`   | Adaptor, links to MAPK cascade |

---

## Basic Queries

### 1.1 Direct interaction partners

Who does the insulin receptor directly interact with?

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

### 1.2 What does IRS1 interact with?

```cypher
MATCH (p:Protein {entryNameID: 'IRS1_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

### 1.3 Shared partners between INSR and IRS1

```cypher
MATCH (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH]-(common:Protein)-[:INTERACTS_WITH]-(p2:Protein {entryNameID: 'IRS1_HUMAN'})
RETURN common.entryNameID AS sharedPartner
```

### 1.4 Shortest path from INSR to AKT1

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'AKT1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 1.5 All shortest paths from INSR to AKT1

```cypher
MATCH path = allShortestPaths(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'AKT1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
ORDER BY hops
```

### 1.6 Most connected proteins in the insulin subgraph

Which proteins are the most central within 3 hops of INSR?

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..3]-(partner:Protein)
WITH partner, count(*) AS connections
RETURN partner.entryNameID, connections
ORDER BY connections DESC
LIMIT 10
```

> Hub proteins like PIK3R1, AKT1, and GRB2 naturally emerge here — they are genuine pathway members but also interact with many other proteins.

### 1.7 Explore the insulin local subgraph

Visualize the interaction neighborhood up to 3 hops from INSR — best run in Neo4j Browser for graph visualization.

```cypher
MATCH path = (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..3]-(partner:Protein)
RETURN path
```

---

## Disease Queries

### 2.1 Diseases associated with the insulin receptor

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:ASSOCIATED_WITH]->(d:Disease)
RETURN d.name, d.id
```

### 2.2 All proteins linked to diabetes

```cypher
MATCH (p:Protein)-[:ASSOCIATED_WITH]->(d:Disease)
WHERE d.name CONTAINS 'diabetes'
RETURN p.entryNameID, p.proteinName, d.name
```

### 2.3 Shortest path between INSR (diabetes) and BRCA1 (cancer)

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..6]-(p2:Protein {entryNameID: 'BRCA1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 2.4 Proteins shared between diabetes and cancer

```cypher
MATCH (d1:Disease)<-[:ASSOCIATED_WITH]-(p:Protein)-[:ASSOCIATED_WITH]->(d2:Disease)
WHERE d1.name CONTAINS 'diabetes' AND d2.name CONTAINS 'cancer'
RETURN p.entryNameID, p.proteinName, d1.name AS disease1, d2.name AS disease2
```

---

## GO Queries

### 3.1 Biological processes of INSR

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INVOLVED_IN]->(bp:BiologicalProcess)
RETURN bp.name, bp.id
```

### 3.2 Proteins sharing a biological process with INSR

```cypher
MATCH (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INVOLVED_IN]->(bp:BiologicalProcess)<-[:INVOLVED_IN]-(p2:Protein)
RETURN bp.name, collect(p2.entryNameID) AS proteins
```

### 3.3 Cellular components of insulin pathway proteins

```cypher
MATCH (p:Protein)-[:IS_LOCATED_IN]->(cc:CellularComponent)
WHERE p.entryNameID IN ['INSR_HUMAN', 'IRS1_HUMAN', 'AKT1_HUMAN', 'MTOR_HUMAN']
RETURN p.entryNameID, collect(cc.name) AS locations
```
