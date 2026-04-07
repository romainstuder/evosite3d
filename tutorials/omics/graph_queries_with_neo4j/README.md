# Neo4j & UniProt: Exploring Human Protein Interactions

A hands-on tutorial for building and querying a human protein-protein interaction (PPI) graph database using Neo4j and UniProt data.

Graph databases excel at revealing hidden relationships in complex networks. Neo4j has been used across many fields — most famously in investigative journalism, where it played a central role in analyzing the **Panama Papers**: 11.5 million leaked financial documents exposing offshore tax havens. Journalists used Neo4j to map the connections between shell companies, intermediaries, and individuals that would have been impossible to detect in a traditional database. Read more: [Analyzing the Panama Papers with Neo4j](https://neo4j.com/blog/cypher-and-gql/analyzing-panama-papers-neo4j/)

The same graph thinking applies to biology — proteins interact in networks, diseases propagate through pathways, and hidden connections emerge when you query the right relationships.

---

## Table of Contents

1. [Why a Graph Database?](#1-why-a-graph-database)
2. [Installing Neo4j](#2-installing-neo4j)
3. [Downloading UniProt Data](#3-downloading-uniprot-data)
4. [Data Ingestion](#4-data-ingestion)
5. [Graph Model](#5-graph-model)
6. [Key Concepts Covered](#6-key-concepts-covered)
7. [Example Queries](#7-example-queries)

---

## 1. Why a Graph Database?

Protein-protein interaction data is inherently a network — proteins are nodes, interactions are edges. While you _could_ store this in a relational database, it quickly becomes painful:

| Question                           | Relational DB               | Neo4j                             |
| ---------------------------------- | --------------------------- | --------------------------------- |
| Who does protein A interact with?  | `JOIN` on interaction table | `MATCH (a)-[:INTERACTS_WITH]-(b)` |
| Find all proteins within 3 hops    | Recursive CTEs, slow        | `[:INTERACTS_WITH*..3]`, native   |
| Shortest path between two proteins | Complex, expensive          | `shortestPath()`, built-in        |
| Add a new relationship type        | Schema migration            | Just add a new relationship       |

The core issue with relational databases is that **path queries require self-joins** — each additional hop multiplies the query complexity. A 5-hop path query on a 20,000-protein dataset becomes impractical in SQL but is trivial in Cypher.

Graph databases also match how biologists think: _"What is between TLR4 and NF-κB?"_ is naturally a path query, not a table scan.

---

## 2. Installing Neo4j

See [INSTALL.md](INSTALL.md) for full installation instructions (Desktop, Server, Docker) and index creation.

---

## 3. Downloading UniProt Data

We use the **UniProt reviewed (Swiss-Prot)** database filtered for _Homo sapiens_ (taxon 9606).

### Steps

1. Go to [https://www.uniprot.org](https://www.uniprot.org)
2. Search for: `reviewed:true AND organism_id:9606`
3. Click **Download**
4. Choose format: **TSV**
5. Select the following columns (`Entry` is selected by default):

| Column                               | Description                                          |
| ------------------------------------ | ---------------------------------------------------- |
| `Entry`                              | UniProt accession (unique ID, e.g. `P04637`)         |
| `Entry Name`                         | Human-readable ID (e.g. `P53_HUMAN`)                 |
| `Protein names`                      | Full protein name                                    |
| `Gene Names`                         | Associated gene names                                |
| `Interacts with`                     | Interacting protein accessions (semicolon-separated) |
| `Involvement in disease`             | Associated diseases                                  |
| `Gene Ontology (biological process)` | GO biological process terms                          |
| `Gene Ontology (molecular function)` | GO molecular function terms                          |
| `Gene Ontology (cellular component)` | GO cellular component terms                          |

6. Download and save as `uniprotkb_human.tsv.gz`
7. Uncompress the file `gunzip uniprotkb_human.tsv.gz`
8. Check number of lines `wc -l uniprotkb_human.tsv` => 20432 uniprotkb_human.tsv

---

## 4. Data Ingestion

Start neo4j:

```bash
neo4j start
```

or

```bash
neo4j restart
```

and access Neo4j Browser at `http://localhost:7474` — default credentials are `neo4j / neo4j`.

First, we need to create indices, which will produce speed:

```cypher
CREATE INDEX protein_index_name FOR (p:Protein) ON (p.entryID);
```

Then the following query will parse the uniprot file and create Protein nodes, with two
properties: `entryID` (based on `Entry` column) and `entryName` (based on `Entry Name` column).

```cypher
WITH "file:///Users/romainstuder/Github/evosite3d/tutorials/omics/graph_queries_with_neo4j/uniprotkb_human.tsv" AS file_path
LOAD CSV WITH HEADERS FROM file_path AS row
FIELDTERMINATOR '\t'
CREATE (p:Protein {entryID:row.Entry, entryNameID:row.`Entry Name`});
```

Note: if you need to delete the current data, you can run:

```cypher
DROP INDEX protein_index_name IF EXISTS;
MATCH (n)
OPTIONAL MATCH (n)-[r]-()
DELETE n, r;
```

You can load some proteins at random:

```cypher
MATCH (n:Protein) RETURN n LIMIT 25;
```

You can search for a specific node using

```cypher
MATCH (n:Protein {entryNameID:"TLR4_HUMAN"}) RETURN *;
```

Now you can load more columns from the uniprot files using the following script:

```bash
cat script.cql | cypher-shell -a bolt://localhost:7687 -u neo4j
```

The ingestion script handles:

- Creating one `(:Protein)` node per row using `Entry` as the unique ID
- Splitting the `Interacts with` column on `";"` and creating `[:INTERACTS_WITH]` relationships
- Parsing `Involvement in disease` entries (format: `DISEASE: <name> [MIM:<id>]: ...`) to create `(:Disease)` nodes with MIM ID and name
- Parsing Gene Ontology columns to create `(:BiologicalProcess)`, `(:MolecularFunction)`, and `(:CellularComponent)` nodes with GO ID and name
- Skipping null or empty fields

---

## 5. Graph Model

```
(:Protein {
    entryID:      "P04637",        // UniProt accession
    entryNameID:  "P53_HUMAN",     // Entry name
    proteinName:  "Cellular tumor antigen p53",
    geneName:     "TP53"
})

(:Disease {
    id:    "151623",               // MIM ID
    name:  "Li-Fraumeni syndrome"  // Disease name
})

(:BiologicalProcess { id: "GO:0006915", name: "apoptotic process" })
(:MolecularFunction { id: "GO:0003677", name: "DNA binding" })
(:CellularComponent { id: "GO:0005634", name: "nucleus" })
```

### Relationships

```
(:Protein)-[:INTERACTS_WITH]->(:Protein)
(:Protein)-[:ASSOCIATED_WITH]->(:Disease)
(:Protein)-[:INVOLVED_IN]->(:BiologicalProcess)
(:Protein)-[:HAS_FUNCTION]->(:MolecularFunction)
(:Protein)-[:IS_LOCATED_IN]->(:CellularComponent)
```

Relationships are directional as loaded from UniProt, but are queried without direction (undirected) unless otherwise specified.

---

## 6. Key Concepts Covered

| Concept                | Cypher Feature                                   |
| ---------------------- | ------------------------------------------------ |
| Single node lookup     | `MATCH (p:Protein {entryNameID: '...'})`         |
| Direct neighbors       | `-[:INTERACTS_WITH]-`                            |
| Aggregation & ranking  | `count()`, `ORDER BY`, `LIMIT`                   |
| 2-hop patterns         | `-[:INTERACTS_WITH]-(common)-[:INTERACTS_WITH]-` |
| Shortest path          | `shortestPath()` with `*..n` bound               |
| All shortest paths     | `allShortestPaths()` with `*..n` bound           |
| Subgraph visualization | `RETURN path` in Neo4j Browser                   |
| Text filtering         | `CONTAINS`, `IS NOT NULL`                        |
| List operations        | `split()`, `size()`                              |

---

## 7. Example Queries

- [Insulin Signaling Pathway](examples_insulin.md) — INSR, IRS1, AKT1, mTOR, disease & GO queries
- [Innate Immunity Signaling](examples_innate_immunity.md) — TLR4, NF-kB, JUN, IRF3, disease & GO queries

---

## Further Reading

- UniProt documentation: [https://www.uniprot.org/help](https://www.uniprot.org/help)
- Cypher query language: [https://neo4j.com/docs/cypher-manual](https://neo4j.com/docs/cypher-manual)
- KEGG TLR signaling pathway: `hsa04620`
- Reactome insulin signaling: [https://reactome.org](https://reactome.org)
