# Structure Server - MCP Tools Documentation

本ドキュメントは `structure_server.py` に実装されているMCPツール関数について、開発者向けの実装詳細と研究者向けの計算手法を説明します。

---

## 目次

1. [fetch_molecules](#1-fetch_molecules)
2. [inspect_molecules](#2-inspect_molecules)
3. [split_molecules](#3-split_molecules)
4. [clean_protein](#4-clean_protein)
5. [clean_ligand](#5-clean_ligand)
6. [run_antechamber_robust](#6-run_antechamber_robust)
7. [prepare_complex](#7-prepare_complex)
8. [create_mutated_structutre](#8-create_mutated_structutre)

---

## 1. fetch_molecules

**概要**: PDB、AlphaFold、PDB-REDOからタンパク質構造ファイルを取得する

### 開発者向け実装詳細

#### 処理フロー

1. **入力検証**
   - `source` パラメータの検証（'pdb', 'alphafold', 'pdb-redo'）
   - PDB IDの大文字変換

2. **ファイル取得** (httpx非同期クライアント使用)
   - **PDB**: `https://files.rcsb.org/download/{pdb_id}.cif` (優先) → `.pdb` (フォールバック)
   - **AlphaFold**: `https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb`
   - **PDB-REDO**: `https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb`

3. **ファイル保存**: `WORKING_DIR/{pdb_id}.{cif|pdb}`

4. **構造統計の取得** (gemmiライブラリ使用)
   - 原子数のカウント
   - チェーンID（label_asym_id）の抽出
   - gemmiが利用不可の場合は簡易パーサーにフォールバック

#### エラーハンドリング

- HTTP 404: 「Structure not found」メッセージとヒント
- タイムアウト: 30秒のタイムアウト設定
- 接続エラー: ネットワーク状態の確認を促すヒント

### 研究者向け計算手法

#### 使用データソース

| ソース | 説明 | フォーマット |
|--------|------|-------------|
| **RCSB PDB** | Protein Data Bankの公式構造データベース | mmCIF (優先), PDB |
| **AlphaFold DB** | DeepMindのAI予測構造データベース | PDB |
| **PDB-REDO** | 結晶構造の自動再精密化データベース | PDB |

#### 出力ファイル形式

- **mmCIF**: 現代的な構造ファイル形式、大規模複合体に対応
- **PDB**: 従来形式、一部のツールとの互換性のため提供

---

## 2. inspect_molecules

**概要**: mmCIF/PDBファイルの分子構成を解析し、チェーン情報を返す

### 開発者向け実装詳細

#### 処理フロー

1. **ファイル読み込み** (gemmi使用)
   ```python
   # mmCIF
   doc = gemmi.cif.read(path)
   structure = gemmi.make_structure_from_block(doc[0])
   
   # PDB
   structure = gemmi.read_pdb(path)
   structure.setup_entities()
   ```

2. **ヘッダー情報の抽出**
   - PDB ID (`structure.name`)
   - 分解能 (`structure.resolution`)
   - 空間群 (`structure.spacegroup_hm`) → X線回折の判定
   - モデル数 > 1 → NMRの判定

3. **エンティティ情報の抽出**
   - `structure.entities` からポリマー種別（polypeptide, polyribonucleotide等）を取得
   - エンティティ名（分子名）の抽出

4. **チェーン分類** (各subchainに対して)
   ```python
   AMINO_ACIDS = {'ALA', 'ARG', ...}  # 22種
   WATER_NAMES = {'HOH', 'WAT', 'H2O', 'DOD', 'D2O'}
   COMMON_IONS = {'NA', 'CL', 'K', 'MG', 'CA', 'ZN', ...}
   
   # 分類ロジック
   if has_protein_residues: chain_type = "protein"
   elif has_water: chain_type = "water"
   elif has_ion: chain_type = "ion"
   else: chain_type = "ligand"
   ```

5. **配列抽出**: タンパク質チェーンの場合、1文字アミノ酸コードの配列を生成

### 研究者向け計算手法

#### 分子分類基準

| 分類 | 基準 |
|------|------|
| **タンパク質** | 標準20アミノ酸 + セレノシステイン(SEC) + ピロリシン(PYL) を含む |
| **水分子** | HOH, WAT, H2O, DOD, D2O |
| **イオン** | NA, CL, K, MG, CA, ZN, FE, MN, CU, CO, NI, CD, HG |
| **リガンド** | 上記以外の低分子化合物 |

#### 出力情報

- **label_asym_id**: mmCIF形式のユニークチェーン識別子（26文字制限なし）
- **auth_asym_id**: 著者指定のチェーンID（従来のA, B, C等）
- **entity情報**: PDBヘッダーからの分子名・記述

---

## 3. split_molecules

**概要**: 構造ファイルを分子種別ごとに個別のPDBファイルに分割する

### 開発者向け実装詳細

#### 処理フロー

1. **構造解析**: `inspect_molecules()` を内部呼び出しして構造情報を取得

2. **チェーン選択**
   - `select_chains` 指定時: 指定チェーンのみ抽出
   - 指定なし: protein + ligand + ion チェーンを抽出（水は除外）

3. **サブチェーン抽出** (gemmi使用)
   ```python
   for subchain in model.subchains():
       chain_id = subchain.subchain_id()  # label_asym_id
       
       # 新規構造を構築
       new_structure = gemmi.Structure()
       new_chain = gemmi.Chain(chain_id)
       
       for residue in subchain:
           # 残基・原子をコピー
           # altloc処理: '\x00', '', 'A', ' ' のみ採用
   ```

4. **ファイル命名規則**
   - `protein_1.pdb`, `protein_2.pdb`, ...
   - `ligand_1.pdb`, `ligand_2.pdb`, ...
   - `ion_1.pdb`, `ion_2.pdb`, ...
   - `water_1.pdb`, `water_2.pdb`, ...

5. **メタデータ保存**: `split_metadata.json` に分割情報を記録

#### Altloc（代替配座）処理

- 結晶構造では同一原子の複数配座が存在する場合がある
- 本実装ではaltloc 'A' または空（デフォルト配座）のみを採用
- 重複原子名は最初の出現のみ保持

### 研究者向け計算手法

#### 出力ファイル形式

すべてPDB形式で出力（mmCIF入力の場合も変換）

#### チェーン識別

mmCIF形式のlabel_asym_idを使用することで、26チェーン以上の大規模複合体にも対応

---

## 4. clean_protein

**概要**: タンパク質構造をMDシミュレーション用に前処理する

### 開発者向け実装詳細

#### 処理フロー（9ステップ）

1. **構造読み込み** (OpenMM PDBFixer)
   ```python
   fixer = PDBFixer(filename=pdb_file)
   ```

2. **欠損残基の検出**
   ```python
   fixer.findMissingResidues()
   # ignore_terminal_missing_residues=True の場合、末端欠損は除外
   ```

3. **末端キャッピング** (オプション)
   ```python
   if cap_termini:
       fixer.missingResidues[chain_idx, 0] = ['ACE']  # N末端
       fixer.missingResidues[chain_idx, length] = ['NME']  # C末端
   ```

4. **非標準残基の処理**
   ```python
   fixer.findNonstandardResidues()
   fixer.replaceNonstandardResidues()  # 標準残基に置換
   ```

5. **ヘテロ原子の除去**
   ```python
   fixer.removeHeterogens(keepWater=keep_water)
   ```

6. **欠損原子の追加**
   ```python
   fixer.findMissingAtoms()
   fixer.addMissingAtoms()  # 重原子・末端原子を追加
   ```

7. **ジスルフィド結合の検出・処理**
   ```python
   fixer.topology.createDisulfideBonds(fixer.positions)
   # CYS → CYX へリネーム（Amber互換性）
   ```

8. **水素原子の付加**
   ```python
   fixer.addMissingHydrogens(pH=ph)  # 指定pHでプロトネーション
   ```

9. **pdb4amberによる命名規則変換**
   ```bash
   pdb4amber -i input.clean.pdb -o output.amber.pdb
   ```

#### 出力ファイル

- `{stem}.clean.pdb`: PDBFixer出力（中間ファイル）
- `{stem}.amber.pdb`: pdb4amber変換後（最終出力）

### 研究者向け計算手法

#### 使用ソフトウェア

| ツール | 目的 | 計算手法 |
|--------|------|----------|
| **OpenMM PDBFixer** | 構造修復・水素付加 | 力場パラメータによる原子配置 |
| **pdb4amber** | Amber命名規則変換 | 残基名・原子名の標準化 |

#### プロトネーション状態

pHに基づいて以下の残基のプロトネーション状態を決定：

- **HIS**: HID (Nδ-H), HIE (Nε-H), HIP (両方プロトン化)
- **ASP/GLU**: 脱プロトン化（pH > pKa）または中性
- **LYS/ARG**: プロトン化（正電荷）

#### ジスルフィド結合

- SG原子間距離が閾値以下のCYS残基を検出
- CYS → CYX へリネームしてAmber力場での正しい処理を保証

---

## 5. clean_ligand

**概要**: リガンド構造を結合次数・プロトネーション状態を正確に設定して前処理する

### 開発者向け実装詳細

#### 処理フロー（8ステップ）

1. **SMILES取得**
   ```python
   # 優先順位: ユーザー指定 > PDB CCD API > 内蔵辞書
   smiles = _get_ligand_smiles(ligand_id, user_smiles, fetch_from_ccd)
   ```

2. **pH依存プロトネーション** (Dimorphite-DL使用)
   ```python
   protonated_smiles, charge = _apply_ph_protonation(smiles, target_ph)
   ```

3. **PDB読み込み** (RDKit)
   ```python
   pdb_mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
   ```

4. **結合次数割り当て** (テンプレートマッチング)
   ```python
   mol_with_bonds = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
   ```

5. **水素原子追加**
   ```python
   mol_with_h = Chem.AddHs(mol_with_bonds, addCoords=True)
   ```

6. **構造最適化** (オプション、MMFF94力場)
   ```python
   AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
   ```

7. **正味電荷の決定**
   - Dimorphite-DL計算結果を優先
   - 手動指定による上書きも可能

8. **SDF出力**
   ```python
   writer = Chem.SDWriter(output_path)
   writer.write(mol_with_h)
   ```

#### SMILES取得元

1. **PDB Chemical Component Dictionary (CCD) API**
   ```
   https://files.rcsb.org/ligands/download/{ligand_id}_ideal.sdf
   → SMILES抽出
   ```

2. **内蔵辞書** (`KNOWN_LIGAND_SMILES`)
   - ATP, ADP, NAD, FAD, SAH等の一般的リガンド

### 研究者向け計算手法

#### 結合次数問題の解決

PDBファイルは結合次数情報を持たないため、以下の手法で正確な結合次数を割り当て：

1. **SMILESテンプレートマッチング**
   - RDKitの `AssignBondOrdersFromTemplate` を使用
   - 3D座標を保持したまま結合次数を転写

2. **PDB CCD (Chemical Component Dictionary)**
   - wwPDBが管理する化学成分の標準定義
   - 正確なSMILES、立体化学情報を提供

#### プロトネーション状態計算

**Dimorphite-DL** を使用したpH依存プロトネーション：

- 各イオン化可能基のpKa予測に基づく
- 生理的pH (7.4) での主要プロトネーション状態を生成
- 複数のプロトネーション状態がある場合は最も確率の高いものを採用

#### 構造最適化

**MMFF94 (Merck Molecular Force Field)**
- 低分子に最適化された力場
- 結合長・結合角・ねじれ角・非結合相互作用を最適化
- 最大200イテレーションまたは収束まで

---

## 6. run_antechamber_robust

**概要**: Antechamberを使用してリガンドのGAFF/GAFF2力場パラメータを生成する

### 開発者向け実装詳細

#### 処理フロー

1. **電荷自動推定** (未指定時)
   ```python
   charge_result = estimate_net_charge(ligand_file)
   net_charge = charge_result["estimated_charge_at_ph"]
   ```

2. **Antechamber実行**
   ```bash
   antechamber -i ligand.sdf -fi mdl \
               -o ligand.gaff.mol2 -fo mol2 \
               -c bcc -nc {charge} \
               -at gaff2 -rn LIG -pf y -j 5
   ```

3. **リトライ機構**
   - 失敗時に電荷を±1調整して再試行
   - 接続性問題の場合、OpenBabelで結合を再構築

4. **parmchk2実行**
   ```bash
   parmchk2 -i ligand.gaff.mol2 -f mol2 \
            -o ligand.frcmod -s gaff2
   ```

5. **frcmod検証**
   - "ATTN" マークのパラメータをチェック
   - 不足パラメータの警告を生成

6. **部分電荷抽出** (MOL2から)
   ```python
   # @<TRIPOS>ATOM セクションから電荷値を抽出
   charges = [float(parts[8]) for line in atom_section]
   ```

#### 診断情報

- `sqm.out`: 半経験的量子計算の出力
- `sqm.in`: SQM入力ファイル
- `diagnostics/`: 各試行の診断ファイル

### 研究者向け計算手法

#### 使用ソフトウェア (AmberTools)

| ツール | 目的 | 計算手法 |
|--------|------|----------|
| **Antechamber** | 原子タイプ・電荷割り当て | AM1-BCC電荷計算 |
| **parmchk2** | 力場パラメータ補完 | 類似パラメータからの推定 |
| **SQM** | 半経験的QM計算 | AM1法 |

#### AM1-BCC電荷計算

1. **AM1 (Austin Model 1)**
   - 半経験的分子軌道法
   - Mulliken電荷を計算

2. **BCC (Bond Charge Correction)**
   - 結合タイプに基づく電荷補正
   - HF/6-31G* RESP電荷を再現するよう最適化

#### GAFF2力場

**General Amber Force Field 2**
- 医薬品様分子に最適化
- 結合伸縮、角度変角、二面角回転、非結合相互作用のパラメータ
- GAFFの改良版（より正確なパラメータ）

#### frcmodファイル

力場修正ファイル - 標準パラメータにない結合・角度・二面角のパラメータを追加

- "ATTN" マーク: 類似パラメータからの推定値（要注意）
- "same as" マーク: 他のパラメータを流用

---

## 7. prepare_complex

**概要**: タンパク質-リガンド複合体の完全な前処理ワークフローを一括実行

### 開発者向け実装詳細

#### 処理フロー（5ステップ統合）

1. **構造解析**: `inspect_molecules()` で分子構成を把握

2. **構造分割**: `split_molecules()` でチェーン別ファイル生成

3. **タンパク質処理** (各タンパク質チェーンに対して)
   ```python
   for protein_file in split_result['protein_files']:
       clean_result = clean_protein(protein_file, ph=ph, cap_termini=cap_termini)
   ```

4. **リガンド処理** (各リガンドチェーンに対して)
   ```python
   for ligand_file in split_result['ligand_files']:
       clean_result = clean_ligand(ligand_file, ligand_id, target_ph=ph)
       if run_parameterization:
           param_result = run_antechamber_robust(clean_result['sdf_file'], ...)
   ```

5. **サマリー生成**: `prepare_complex_summary.json`

#### 出力構造

```
output_dir/
├── protein_1.pdb
├── protein_1.clean.pdb
├── protein_1.amber.pdb
├── ligand_1.pdb
├── ligand_1_prepared.sdf
├── ligand_1_prepared.gaff.mol2
├── ligand_1_prepared.frcmod
├── split_metadata.json
└── prepare_complex_summary.json
```

### 研究者向け計算手法

#### ワークフロー概要

```
入力 (mmCIF/PDB)
    ↓
[inspect_molecules] 構造解析
    ↓
[split_molecules] チェーン分割
    ↓
┌─────────────────┬─────────────────┐
↓                 ↓                 ↓
タンパク質      リガンド          イオン
    ↓                 ↓
[clean_protein]   [clean_ligand]
• PDBFixer        • SMILES取得
• 水素付加        • Dimorphite-DL
• pdb4amber       • RDKit最適化
    ↓                 ↓
.amber.pdb      .sdf
                      ↓
              [run_antechamber_robust]
              • AM1-BCC電荷
              • GAFF2パラメータ
                      ↓
              .mol2 + .frcmod
```

#### 対応入力形式

- **実験構造**: PDB/mmCIF (X線、NMR、クライオEM)
- **計算構造**: Boltz-2、AlphaFold等のAI予測構造

---

## 8. create_mutated_structutre

**概要**: タンパク質構造に点変異を導入する

### 開発者向け実装詳細

#### 処理フロー

1. **配列抽出** (`pdb_to_sequence`)
   ```python
   # CA原子から1文字アミノ酸配列を構築
   for res in structure.residues():
       if 'CA' in res.atoms():
           sequence.append(AA_CODE[res.name])
   ```

2. **変異指定の解析** (`create_mutation_dict`)
   ```python
   # "10,25,100" + "A,G,W" → {10: 'A', 25: 'G', 100: 'W'}
   mutation_dict = create_mutation_dict(mutation_indices, mutation_residues)
   ```

3. **変異配列の生成**
   ```python
   mutated_sequence = sequence.copy()
   mutated_sequence[idx - 1] = new_aa
   ```

4. **FASPR側鎖パッキング** (`generate_structure`)
   - 変異残基の側鎖配座を予測
   - 周囲残基との立体衝突を最小化

5. **PDBFixer仕上げ**
   ```python
   fixer = PDBFixer(mutated_pdb)
   PDBFile.writeFile(fixer.topology, fixer.positions, output)
   ```

### 研究者向け計算手法

#### 使用ソフトウェア

| ツール | 目的 | 計算手法 |
|--------|------|----------|
| **FASPR** | 側鎖パッキング | ロータマーライブラリベースの最適化 |
| **PDBFixer** | 構造クリーニング | トポロジー修正 |

#### FASPR (Fast Side-chain Packing using Rotamers)

- **ロータマーライブラリ**: 各アミノ酸の統計的に好ましい側鎖配座
- **エネルギー最小化**: 
  - 立体衝突（van der Waals反発）
  - 静電相互作用
  - 水素結合
- **計算効率**: グラフ理論に基づく高速アルゴリズム

#### 変異表記法

- **標準表記**: `{元の残基}{位置番号}{変異後残基}` (例: A100G)
- **位置番号**: 1-based (PDBの残基番号ではなく配列位置)

---

## 依存ライブラリ

### Python パッケージ

| パッケージ | バージョン | 用途 |
|------------|-----------|------|
| gemmi | >= 0.6.0 | mmCIF/PDB解析 |
| openmm | >= 8.0 | PDBFixer |
| rdkit | >= 2023.09 | 化学構造処理 |
| httpx | >= 0.24 | 非同期HTTP |
| dimorphite-dl | >= 1.3 | プロトネーション |

### 外部コマンド (AmberTools)

| コマンド | 用途 |
|----------|------|
| pdb4amber | PDB命名規則変換 |
| antechamber | リガンドパラメータ化 |
| parmchk2 | 力場パラメータ補完 |
| sqm | 半経験的QM計算 |

### オプション

| ツール | 用途 |
|--------|------|
| FASPR | 側鎖パッキング |
| OpenBabel | 分子ファイル変換 |

---

## エラーハンドリング

すべてのMCPツールは以下の構造で結果を返します：

```python
{
    "success": bool,      # 処理成功/失敗
    "errors": list[str],  # エラーメッセージ（失敗時）
    "warnings": list[str], # 警告メッセージ
    "output_file": str,   # 出力ファイルパス
    # ... ツール固有のフィールド
}
```

### エラーメッセージの形式

- **エラー本文**: `{ErrorType}: {message}`
- **ヒント**: `Hint: {解決策の提案}`

これにより、LLMエージェントが自動的にエラーを理解し、適切な対処を行えます。

---

## 参考文献

### ソフトウェア

1. **PDBFixer**: Eastman et al., "OpenMM 7: Rapid development of high performance algorithms for molecular dynamics", PLoS Comput Biol (2017)

2. **RDKit**: Open-source cheminformatics. https://www.rdkit.org

3. **gemmi**: Wojdyr, "GEMMI: A library for structural biology", J. Open Source Softw. (2022)

4. **AmberTools**: Case et al., "Amber 2023", University of California, San Francisco

5. **Dimorphite-DL**: Ropp et al., "Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules", J. Cheminform. (2019)

6. **FASPR**: Huang et al., "FASPR: an open-source tool for fast and accurate protein side-chain packing", Bioinformatics (2020)

### 力場・電荷モデル

7. **GAFF2**: He et al., "A fast and high-quality charge model for the next generation general AMBER force field", J. Chem. Phys. (2020)

8. **AM1-BCC**: Jakalian et al., "Fast, efficient generation of high-quality atomic charges. AM1-BCC model: II. Parameterization and validation", J. Comput. Chem. (2002)

### データベース

9. **RCSB PDB**: Berman et al., "The Protein Data Bank", Nucleic Acids Res. (2000)

10. **AlphaFold Protein Structure Database**: Varadi et al., "AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models", Nucleic Acids Res. (2022)

11. **PDB Chemical Component Dictionary**: Westbrook et al., "The Chemical Component Dictionary: complete descriptions of constituent molecules in experimentally determined 3D macromolecular structures", Bioinformatics (2015)

