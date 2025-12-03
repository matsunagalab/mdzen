# Amber Prep Server - 設計仕様書

mmCIF/PDB構造ファイルからAmber/OpenMM用のMDシミュレーション入力ファイルを生成するMCPサーバー。

## 対応入力ソース

- **PDBデータベース**: mmCIF (.cif) または PDB (.pdb) 形式
- **Boltz-2予測出力**: mmCIF形式
- **AlphaFold予測出力**: mmCIF/PDB形式
- **その他の構造予測/実験構造ファイル**

## 概要

### 推奨ワークフロー（SMILESテンプレート法）

```
Structure File (mmCIF/PDB)
       │
       ▼
┌─────────────────────┐
│   parse_structure   │ ← gemmiで解析、チェーン選択、リガンド抽出
└─────────────────────┘
       │
       ├─────────────────────────────────────────────────┐
       ▼                                                 ▼
┌──────────────┐                        ┌───────────────────────────────┐
│prepare_protein│                        │ prepare_ligand_for_amber     │ ★推奨
│  (pdb4amber)  │                        │ (SMILESテンプレート+MMFF94)   │
└──────────────┘                        └───────────────────────────────┘
                                                         │
                                         ┌───────────────┼───────────────┐
                                         ▼               ▼               ▼
                                    CCD API取得    辞書参照     手動指定
                                    (RCSB)        (KNOWN_)    (引数)
                                         └───────────────┼───────────────┘
                                                         ▼
                                              ┌─────────────────────┐
                                              │ AssignBondOrders    │ ← 結合次数確定
                                              │ FromTemplate        │
                                              └─────────────────────┘
                                                         │
                                                         ▼
                                              ┌─────────────────────┐
                                              │ AddHs + MMFF94最適化 │
                                              └─────────────────────┘
                                                         │
                                                         ▼
                                              ┌─────────────────────┐
                                              │ SDF出力 (結合情報保持)│
                                              └─────────────────────┘
                                                         │
                                                         ▼
                                        ┌──────────────────────────┐
                                        │ run_antechamber_robust   │ ← -fi sdf
                                        │ (GAFF2 + AM1-BCC)        │
                                        └──────────────────────────┘
                                                         │
                                                         ▼
                                        ┌──────────────────────────┐
                                        │ validate_frcmod          │
                                        └──────────────────────────┘
                                                         │
                                                         ▼
                                        ┌──────────────────────────┐
                                        │ build_multi_ligand_system│
                                        └──────────────────────────┘
```

### なぜSMILESテンプレート法か？

PDBファイルには**結合次数情報が含まれていません**。従来のOpenBabel/antechamberによる幾何学的な結合推定は、以下の問題を引き起こします：

| 問題 | 症状 | 
|-----|------|
| 芳香環の誤認識 | ベンゼン環がシクロヘキサンとして扱われる |
| 二重結合の消失 | ケトン(C=O)がアルコール(C-O)になる |
| antechamberエラー | "more than one unit" エラー |
| 電荷計算失敗 | 不正な原子タイプ割り当て |

**SMILESテンプレート法**は、PDB Chemical Component Dictionary (CCD) から「正しい化学構造」を取得し、座標にマッピングすることでこれらの問題を**100%排除**します。

### レガシーワークフロー（参考）

```
┌─────────────────────┐     ┌─────────────────┐
│ estimate_net_charge │ ──► │prepare_ligand_h │ ← OpenBabel
│ (RDKit自動推定)     │     │ (水素付加)       │
└─────────────────────┘     └─────────────────┘
```

`prepare_ligand_hydrogens`は引き続き利用可能ですが、複雑なリガンドでは`prepare_ligand_for_amber`の使用を推奨します。

## MCPツール一覧

### 1. `parse_structure` (推奨)
**汎用構造パーサー**: mmCIF/PDBファイルを解析し、チェーン選択とリガンド抽出を行う。

**入力:**
- `structure_file`: mmCIF (.cif) または PDB (.pdb) ファイルパス
- `output_dir`: 出力ディレクトリ（省略時は自動生成）
- `select_chains`: 抽出するチェーンIDのリスト（例: `["A"]`）。省略時は全タンパク質チェーン
- `include_ligands`: 選択チェーンに結合したリガンドを含める（デフォルト: True）
- `exclude_waters`: 結晶水を除外（デフォルト: True）
- `ligand_distance_cutoff`: リガンドがチェーンに「結合している」とみなす距離閾値（デフォルト: 5.0 Å）

**出力:**
- `protein_pdb`: タンパク質PDBパス
- `ligand_files`: リガンドPDBパスのリスト
- `all_chains`: 全チェーン情報
- `selected_chains`: 選択されたチェーン
- `extracted_ligands`: 抽出されたリガンド情報

**ファイル命名規則:**
```
output/amber_prep/{job_id}/
├── protein.pdb
├── ligand_AP5_chainA.pdb     # チェーンAに結合したリガンド
├── parse_metadata.json       # メタデータ
```

**使用例:**
```python
# PDBからダウンロードした1AKEのチェーンAとそのリガンドを抽出
result = parse_structure(
    "1AKE.cif",
    select_chains=["A"],
    include_ligands=True,
    exclude_waters=True
)

# 全チェーンを抽出（デフォルト）
result = parse_structure("structure.pdb")

# 特定チェーンのみ、リガンドなし
result = parse_structure(
    "complex.cif",
    select_chains=["A", "B"],
    include_ligands=False
)
```

### 2. `prepare_ligand_for_amber` ★推奨

**SMILESテンプレート法 + Dimorphite-DLプロトン化**によるリガンド準備。PDB CCDから正しい化学構造を取得し、pH依存のプロトン化状態を決定する最も堅牢な方法。

**ワークフロー:**
1. SMILES取得（CCD API > 辞書 > 手動指定）
2. **Dimorphite-DLでpH 7.4プロトン化**（中性SMILES → イオン化状態）
3. テンプレートマッチングで結合次数確定
4. 水素付加 + MMFF94最適化
5. SDFで出力（結合次数保持）
6. **プロトン化分子から正確な総電荷を計算**

**入力:**
- `ligand_pdb`: リガンドPDBファイル（parse_structureから）
- `ligand_id`: 3文字のリガンド残基名（例: "ATP", "SAH"）
- `smiles`: 手動指定SMILES（最優先、省略可）
- `output_dir`: 出力ディレクトリ（省略時は入力ファイルと同じ）
- `optimize`: MMFF94最適化を実行（デフォルト: True）
- `max_opt_iters`: 最適化の最大イテレーション（デフォルト: 200）
- `fetch_smiles`: CCD APIからSMILESを取得（デフォルト: True）
- `target_ph`: プロトン化のターゲットpH（デフォルト: 7.4）★新規
- `manual_charge`: 手動指定の総電荷（Dimorphite-DLが失敗した場合のオーバーライド）★新規

**出力:**
- `sdf_file`: 準備済みSDFファイル（antechamberに直接渡す）
- `net_charge`: 計算された総電荷（pH依存、常に0とは限らない）
- `charge_source`: 電荷の計算元（"dimorphite", "manual"）★新規
- `mol_formal_charge`: RDKit分子の形式電荷★新規
- `smiles_used`: 使用されたSMILES（プロトン化後）
- `smiles_original`: 元のSMILES（プロトン化前）★新規
- `smiles_source`: SMILESの取得元（"user", "ccd", "dictionary"）
- `target_ph`: 使用されたpH★新規
- `optimized`: 最適化実行の有無
- `optimization_converged`: 最適化が収束したか

**SMILESの優先順位:**
1. `smiles`引数（ユーザー指定）
2. PDB CCD API (`https://data.rcsb.org/rest/v1/core/chemcomp/{id}`)
3. `KNOWN_LIGAND_SMILES`辞書（フォールバック）

**ファイル命名規則:**
```
ligand_AP5_chainA.pdb → ligand_AP5_chainA_prepared.sdf
```

**使用例:**
```python
# CCD APIから自動取得 + Dimorphite-DLでpH 7.4プロトン化
result = prepare_ligand_for_amber(
    "ligand_ATP_chainA.pdb",
    "ATP",
    target_ph=7.4,  # 生理的pH（Dimorphite-DLで処理）
    optimize=True
)
# result['net_charge'] は ATP の場合 -4 程度（リン酸基が脱プロトン化）

# 手動でSMILESと電荷を指定（複雑なリガンドの場合）
result = prepare_ligand_for_amber(
    "ligand_CUSTOM_chainA.pdb",
    "CUS",
    smiles="c1ccc(cc1)C(=O)O",  # ベンゾ酸の例
    manual_charge=-1,  # カルボン酸が脱プロトン化
    optimize=True
)

# antechamberに渡す（電荷が自動的に正しく設定される）
antechamber_result = run_antechamber_robust(
    result['sdf_file'],  # SDFを使用
    net_charge=result['net_charge']  # Dimorphite-DLが計算した正確な電荷
)
```

**Note:** CCD APIのSMILESは通常中性形式で、`GetFormalCharge()`は0を返します。Dimorphite-DLがpH依存のイオン化を適用することで、正確な総電荷が得られます。

### 3. `prepare_ligand_hydrogens`（レガシー）

OpenBabelを使用してリガンドに水素を付加。**複雑なリガンドでは`prepare_ligand_for_amber`を推奨。**

**入力:**
- `ligand_file`: リガンド構造ファイル（PDB/MOL2/SDF）
- `ph`: プロトン化のpH（デフォルト: 7.4）
- `output_format`: 出力形式（デフォルト: mol2）

**出力:**
- `output_file`: 水素付加後のファイル
- `num_hydrogens_added`: 追加された水素数

**ファイル命名規則:**
```
ligand_SAH_chainC.pdb → ligand_SAH_chainC_H.mol2
```

### 4. `estimate_net_charge`
RDKitを使用してリガンドの総電荷を推定。

**入力:**
- `ligand_file`: リガンド構造ファイル
- `ph`: ターゲットpH（デフォルト: 7.4）

**出力:**
- `formal_charge`: 形式電荷
- `estimated_charge_at_ph`: pH考慮した推定電荷
- `ionizable_groups`: イオン化基のリスト
- `confidence`: 推定の信頼度（high/medium/low）

**検出するイオン化基:**
| 基 | pKa範囲 | pH 7.4での電荷 |
|---|---------|---------------|
| カルボン酸 | 3-5 | -1 |
| 一級アミン | 9-11 | +1 |
| 二級アミン | 9-11 | +1 |
| スルホン酸 | <1 | -1 |
| リン酸 | 2, 7 | -2 |
| フェノール | 9-10 | 0 |

### 5. `prepare_protein_for_amber`
pdb4amberを使用してタンパク質をAmber用に準備。

**入力:**
- `pdb_file`: 入力PDBファイル
- `detect_disulfides`: ジスルフィド結合自動検出（デフォルト: True）

**出力:**
- `output_pdb`: 準備済みPDB
- `disulfide_bonds`: 検出したSS結合
- `histidine_states`: HIS→HID/HIE/HIP変換情報

**ファイル命名規則:**
```
protein.pdb → protein_amber.pdb
```

### 6. `run_antechamber_robust`
エラーハンドリング強化版antechamber。sqm失敗時に電荷を±1調整して自動リトライ。

**入力:**
- `ligand_file`: リガンドファイル（水素付加済み推奨）
- `net_charge`: 総電荷（省略時は自動推定）
- `residue_name`: 残基名（デフォルト: "LIG"）
- `charge_method`: 電荷計算法（デフォルト: "bcc" = AM1-BCC）
- `atom_type`: 原子タイプ（デフォルト: "gaff2"）
- `max_retries`: リトライ回数（デフォルト: 2）

**出力:**
- `mol2`: GAFF2パラメータ付きMOL2
- `frcmod`: Force modification file
- `charge_used`: 使用した総電荷
- `sqm_diagnostics`: sqm実行診断情報

**ファイル命名規則:**
```
ligand_SAH_chainC_H.mol2 → ligand_SAH_chainC_H.gaff.mol2
                        → ligand_SAH_chainC_H.frcmod
```

**sqmエラー診断:**
| エラーパターン | 診断 | 推奨対処 |
|--------------|------|---------|
| `number of electrons is odd` | 総電荷が不正 | 電荷を±1調整（自動リトライ） |
| `No convergence in SCF` | SCF収束失敗 | 構造最適化後に再試行 |
| `Cannot properly run sqm` | sqm実行失敗 | AmberTools環境確認 |
| `more than one unit` | 結合情報不正 | OpenBabelで結合修復（自動） |

**結合修復機能:**
AP5やATPのような複雑な分子では、PDB/mmCIFから抽出した際に結合情報が不完全になることがあります。
`run_antechamber_robust`は以下の自動修復を行います：
1. `-j 5` オプションで距離ベースの結合追加を試行
2. 失敗した場合、OpenBabelの`--connect`オプションで結合を再構築
3. 修復後のファイルで再試行

### 7. `validate_frcmod`
frcmodファイルの品質検証。

**入力:**
- `frcmod_file`: .frcmodファイルパス

**出力:**
- `valid`: 検証結果
- `attn_count`: "ATTN, need revision"の検出数
- `warnings`: 警告リスト
- `recommendations`: 対処推奨事項

### 8. `build_multi_ligand_system`
複数リガンドを含むMDシステムをtleapで構築。

**入力:**
- `protein_pdb`: 準備済みタンパク質PDB
- `ligands`: リガンド情報のリスト
  ```python
  [
      {'mol2': 'path/to/lig1.gaff.mol2', 'frcmod': 'path/to/lig1.frcmod', 'residue_name': 'SAH'},
      {'mol2': 'path/to/lig2.gaff.mol2', 'frcmod': 'path/to/lig2.frcmod', 'residue_name': 'LIG'},
  ]
  ```
- `water_model`: 水モデル（デフォルト: "tip3p"）
- `box_padding`: ボックスパディング（デフォルト: 12.0 Å）
- `box_type`: ボックス形状（デフォルト: "box"=直方体、"oct"=truncated octahedron）
- `neutralize`: イオンで中和（デフォルト: True）
- `salt_conc`: 塩濃度（デフォルト: 0.15 M）

**出力:**
- `parm7`: Amberトポロジーファイル
- `rst7`: Amber座標ファイル
- `complex_pdb`: 複合体PDB
- `num_atoms`, `num_residues`: 原子・残基数
- `ligand_names`: 含まれるリガンド名

### 9. `build_complex_system`
単一リガンドのMDシステム構築（レガシー互換）。

### 10. `boltz2_to_amber_complete`
Boltz-2 mmCIFからMD入力ファイルまでの完全ワークフロー（単一リガンド用）。

---

## 典型的なワークフロー

### ワークフロー1: PDBデータベースからの構造準備（推奨）

PDBからダウンロードした構造（例: 1AKE）から特定チェーンとそのリガンドを抽出してMD準備。
**SMILESテンプレート法**を使用して堅牢なリガンド準備を行う。

```python
from pathlib import Path
from servers.amber_prep_server import (
    parse_structure,
    prepare_ligand_for_amber,  # ★推奨
    run_antechamber_robust,
    prepare_protein_for_amber,
    build_multi_ligand_system
)

# Step 1: mmCIF/PDB解析（チェーン選択付き）
result = parse_structure(
    "1AKE.cif",
    select_chains=["A"],        # チェーンAのみ抽出
    include_ligands=True,       # 結合リガンドを含む
    exclude_waters=True,        # 結晶水を除外
    ligand_distance_cutoff=5.0  # 5Å以内のリガンドを「結合」とみなす
)
job_dir = result['output_dir']

# Step 2: タンパク質準備
protein_result = prepare_protein_for_amber(
    result['protein_pdb'],
    output_dir=job_dir
)

# Step 3-5: 各リガンドのパラメータ化（SMILESテンプレート法）
ligand_params = []
for lig_file in result['ligand_files']:
    # リガンド名を抽出
    res_name = Path(lig_file).stem.split('_')[1]  # e.g., "AP5"
    
    # ★ SMILESテンプレート法でリガンド準備
    # CCD APIから自動でSMILESを取得、結合次数を確定
    prep_result = prepare_ligand_for_amber(
        ligand_pdb=lig_file,
        ligand_id=res_name,
        optimize=True  # MMFF94で軽く最適化
    )
    
    # antechamber (SDF入力で結合次数を保持)
    ac_result = run_antechamber_robust(
        prep_result['sdf_file'],  # SDFを使用
        net_charge=prep_result['net_charge'],  # SMILESから計算
        residue_name=res_name[:3].upper()
    )
    
    ligand_params.append({
        'mol2': ac_result['mol2'],
        'frcmod': ac_result['frcmod'],
        'residue_name': res_name[:3].upper()
    })

# Step 6: システム構築
system_result = build_multi_ligand_system(
    protein_pdb=protein_result['output_pdb'],
    ligands=ligand_params,
    output_dir=job_dir,
    water_model="tip3p",
    box_padding=10.0,
    box_type="box"  # 直方体
)

print(f"Topology: {system_result['parm7']}")
print(f"Coordinates: {system_result['rst7']}")
```

### ワークフロー2: Boltz-2出力の準備

```python
from pathlib import Path
from servers.amber_prep_server import (
    parse_structure,
    prepare_ligand_for_amber,
    run_antechamber_robust,
    prepare_protein_for_amber,
    build_multi_ligand_system
)

# Step 1: mmCIF解析（全チェーン抽出）
result = parse_structure("boltz_output.cif")
job_dir = result['output_dir']

# Step 2: タンパク質準備
protein_result = prepare_protein_for_amber(
    result['protein_pdb'],
    output_dir=job_dir
)

# Step 3-5: 各リガンドのパラメータ化
ligand_params = []
for lig_file in result['ligand_files']:
    res_name = Path(lig_file).stem.split('_')[1]
    
    # SMILESテンプレート法
    prep_result = prepare_ligand_for_amber(
        ligand_pdb=lig_file,
        ligand_id=res_name,
        optimize=True
    )
    
    ac_result = run_antechamber_robust(
        prep_result['sdf_file'],
        net_charge=prep_result['net_charge'],
        residue_name=res_name[:3].upper()
    )
    
    ligand_params.append({
        'mol2': ac_result['mol2'],
        'frcmod': ac_result['frcmod'],
        'residue_name': res_name[:3].upper()
    })

# Step 6: システム構築
system_result = build_multi_ligand_system(
    protein_pdb=protein_result['output_pdb'],
    ligands=ligand_params,
    output_dir=job_dir,
    water_model="tip3p",
    box_padding=10.0,
    box_type="box"
)

print(f"Topology: {system_result['parm7']}")
print(f"Coordinates: {system_result['rst7']}")
```

### ワークフロー3: 手動SMILESを指定する場合

CCD APIにないカスタムリガンドの場合：

```python
# 手動でSMILESを指定
prep_result = prepare_ligand_for_amber(
    ligand_pdb="custom_ligand.pdb",
    ligand_id="CUS",
    smiles="c1ccc(cc1)C(=O)Nc2ccc(cc2)O",  # 手動指定
    optimize=True
)
```

### ワークフロー4: OpenMMでのシミュレーション実行

```python
import openmm as mm
from openmm import app, unit
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation

# システム読み込み
prmtop = AmberPrmtopFile(system_result['parm7'])
inpcrd = AmberInpcrdFile(system_result['rst7'])

# システム作成
system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=10*unit.angstrom,
    constraints=app.HBonds
)

# NPTアンサンブル
system.addForce(mm.MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin))

# インテグレータ
integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin, 
    1/unit.picosecond, 
    2*unit.femtoseconds
)

# シミュレーション
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# エネルギー最小化
simulation.minimizeEnergy()

# 平衡化・プロダクション
simulation.step(10000)  # 20 ps
```

---

## 出力ディレクトリ構造

```
output/amber_prep/{job_id}/
├── protein.pdb                      # 元のタンパク質
├── protein_amber.pdb                # pdb4amber処理後
├── protein_pdb4amber.log            # pdb4amberログ
│
├── ligand_SAH_chainC.pdb            # 元のリガンド
├── ligand_SAH_chainC_H.mol2         # 水素付加後
├── ligand_SAH_chainC_H.gaff.mol2    # GAFF2パラメータ付き
├── ligand_SAH_chainC_H.frcmod       # 力場修正ファイル
│
├── ligand_LIG1_chainE.pdb
├── ligand_LIG1_chainE_H.mol2
├── ligand_LIG1_chainE_H.gaff.mol2
├── ligand_LIG1_chainE_H.frcmod
│
├── complex.parm7                    # Amberトポロジー
├── complex.rst7                     # Amber座標
├── complex.pdb                      # 複合体PDB（溶媒込み）
├── leap.in                          # tleapスクリプト
├── leap.log                         # tleapログ
│
└── diagnostics/
    ├── sqm_attempt1.out             # sqm出力
    ├── charge_estimation.json       # 電荷推定結果
    └── frcmod_validation.json       # frcmod検証結果
```

---

## 既知の制限事項

### 1. 複雑な分子の電荷推定
SAH（S-adenosyl-L-homocysteine）のような複雑な補酵素は、RDKitでのSMARTSパターンマッチングが正確に機能しない場合があります。

**対処法:** `KNOWN_CHARGES`辞書に手動で電荷を指定

```python
KNOWN_CHARGES = {
    "SAH": -1,   # S-adenosyl-L-homocysteine
    "ATP": -4,   # Adenosine triphosphate
    "ADP": -3,   # Adenosine diphosphate
    "NAD": -1,   # NAD+
    "FAD": -2,   # Flavin adenine dinucleotide
}
```

### 2. 芳香環のケクレ化エラー
RDKitがPDBからの分子読み込み時に「Can't kekulize mol」エラーを出す場合があります（特にプリン環など）。

**対処法:** 
- OpenBabelで水素付加後のMOL2を使用
- 電荷推定が失敗しても、antechamberは内部でsqmを使用するため処理可能

### 3. 同名リガンドの重複
Boltz-2出力で同じ残基名のリガンドが複数チェーンに存在する場合。

**対処法:** ファイル名にチェーン名を含める（自動対応済み）
```
ligand_SAH_chainC.pdb
ligand_SAH_chainD.pdb
```

### 4. Alternate Conformations（構造障害）
PDBファイルに複数のコンフォメーション（altloc A/B）がある場合、antechamberが「more than one unit」エラーを出す。

**対処法:** `parse_structure`でリガンド抽出時に自動的にaltloc ''または'A'のみ保持（自動対応済み）

---

## 依存関係

### Python パッケージ
- `fastmcp` - MCPサーバーフレームワーク
- `gemmi` - mmCIF解析
- `rdkit` - 分子操作、電荷推定
- `py3Dmol` - 3D可視化
- `mdtraj` - トラジェクトリ解析
- `openmm` - MDシミュレーション

### 外部ツール（conda経由）
- `antechamber` - リガンドパラメータ化
- `parmchk2` - パラメータチェック
- `pdb4amber` - タンパク質準備
- `tleap` - システム構築
- `obabel` - 水素付加

### インストール
```bash
# conda環境作成
conda create -n mcp-md python=3.11
conda activate mcp-md

# AmberTools & OpenMM
conda install -c conda-forge ambertools openmm rdkit

# Python依存関係
pip install -e .
```

---

## テスト

`notebooks/test_amber_prep.ipynb`で以下をテスト可能：

1. **Test 5**: Boltz-2 mmCIF解析
2. **Test 6**: Antechamber力場生成
3. **Test 7**: tleapシステム構築（複数リガンド）
4. **Test 8**: tleapビルド結果可視化（parm7/rst7→PDB）
5. **Test 9**: OpenMMシミュレーション
6. **Test 10**: トラジェクトリ可視化

---

---

## 開発者向け情報

### コード構造

```
servers/amber_prep_server.py
├── インポート & 定数定義 (L1-70)
│   ├── FastMCP サーバー初期化
│   ├── WORKING_DIR 定義
│   ├── ツールラッパー初期化 (antechamber, parmchk2, etc.)
│   └── KNOWN_LIGAND_SMILES 辞書  # よく使うリガンドのSMILES
│
├── ヘルパー関数 (L72-340)
│   ├── generate_job_id()         # ジョブID生成 (common.utils から)
│   │
│   ├── # SMILES テンプレート関連 (新規)
│   ├── _fetch_smiles_from_ccd()  # PDB CCDからSMILES取得
│   ├── _get_ligand_smiles()      # SMILES取得（優先順位付きフォールバック）
│   ├── _assign_bond_orders_from_smiles()  # テンプレートマッチング
│   ├── _optimize_ligand_rdkit()  # MMFF94最適化
│   │
│   ├── # sqm/antechamber関連
│   ├── _parse_sqm_output()       # sqm.out解析
│   ├── _parse_frcmod_warnings()  # frcmod検証
│   ├── _estimate_charge_rdkit()  # RDKit電荷推定
│   └── _estimate_physiological_charge()  # pH考慮電荷計算
│
└── MCPツール (L340-2000)
    ├── parse_structure()
    ├── prepare_ligand_for_amber()  # ★新規：SMILESテンプレート法
    ├── prepare_ligand_hydrogens()  # レガシー
    ├── estimate_net_charge()
    ├── prepare_protein_for_amber()
    ├── run_antechamber_robust()
    ├── validate_frcmod()
    ├── build_complex_system()
    ├── build_multi_ligand_system()
    └── boltz2_to_amber_complete()
```

### ツールラッパーの使い方

外部コマンド実行には`common/base.py`の`BaseToolWrapper`を使用:

```python
from common.base import BaseToolWrapper

# ラッパー初期化（conda環境指定）
antechamber_wrapper = BaseToolWrapper("antechamber", conda_env="mcp-md")

# コマンド実行
result = antechamber_wrapper.run(
    args=['-i', 'input.mol2', '-o', 'output.mol2', '-c', 'bcc'],
    cwd=working_dir,  # 作業ディレクトリ（省略可）
    timeout=300       # タイムアウト秒（省略可）
)

# 結果
print(result.stdout)
print(result.stderr)
print(result.returncode)
```

### 新しいMCPツールの追加方法

```python
@mcp.tool()
def my_new_tool(
    input_file: str,
    option1: str = "default",
    option2: Optional[int] = None
) -> dict:
    """ツールの説明（MCPクライアントに表示される）。
    
    Args:
        input_file: 入力ファイルパス
        option1: オプション1の説明
        option2: オプション2の説明（省略可）
    
    Returns:
        Dict with results
    """
    logger.info(f"Running my_new_tool: {input_file}")
    
    # 入力検証
    input_path = Path(input_file).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input not found: {input_file}")
    
    # 出力ディレクトリ設定
    output_dir = input_path.parent
    ensure_directory(output_dir)
    
    # 処理実行
    try:
        # ... 実装 ...
        pass
    except Exception as e:
        logger.error(f"my_new_tool failed: {e}")
        raise
    
    # 結果を辞書で返す（JSON serializable）
    return {
        "input_file": str(input_file),
        "output_file": str(output_file),
        "status": "success",
        # ... その他の結果 ...
    }
```

### ヘルパー関数の拡張

#### 新しいイオン化基の追加

`_estimate_charge_rdkit()`にSMARTSパターンを追加:

```python
def _estimate_charge_rdkit(mol) -> Dict[str, Any]:
    # ... 既存コード ...
    
    # 新しいイオン化基を追加
    # 例: イミダゾール（ヒスチジン側鎖）
    imidazole_pattern = Chem.MolFromSmarts("[nR1]1[cR1][nR1][cR1][cR1]1")
    if imidazole_pattern and mol.HasSubstructMatch(imidazole_pattern):
        matches = mol.GetSubstructMatches(imidazole_pattern)
        result["ionizable_groups"].append({
            "type": "imidazole",
            "count": len(matches),
            "typical_charge": 0,  # pH 7.4ではほぼ中性
            "pka_range": "6.0"
        })
    
    return result
```

#### sqmエラーパターンの追加

`_parse_sqm_output()`に新しいエラーパターンを追加:

```python
def _parse_sqm_output(sqm_out_path: Path) -> Dict[str, Any]:
    # ... 既存コード ...
    
    # 新しいエラーパターン
    if "memory allocation failed" in content.lower():
        diagnostics["errors"].append("Memory allocation failed")
        diagnostics["recommendations"].append(
            "Try reducing molecule size or increasing system memory."
        )
    
    return diagnostics
```

### エラーハンドリングパターン

#### 基本パターン

```python
@mcp.tool()
def robust_tool(input_file: str) -> dict:
    logger.info(f"Starting: {input_file}")
    
    # 1. 入力検証
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"File not found: {input_file}")
    
    # 2. 処理実行（try-except）
    try:
        result = some_operation(input_path)
    except SpecificError as e:
        logger.error(f"Specific error: {e}")
        raise RuntimeError(f"Operation failed: {e}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise
    
    # 3. 出力検証
    if not output_path.exists():
        raise RuntimeError(f"Output not created: {output_path}")
    
    return {"status": "success", ...}
```

#### リトライパターン（antechamber参考）

```python
def robust_operation_with_retry(input_file, max_retries=3):
    params_to_try = [default_param, alt_param1, alt_param2]
    last_error = None
    
    for attempt, param in enumerate(params_to_try[:max_retries]):
        logger.info(f"Attempt {attempt + 1}: param = {param}")
        
        try:
            result = operation(input_file, param)
            if verify_success(result):
                logger.info(f"Succeeded with param = {param}")
                return result
        except Exception as e:
            last_error = e
            logger.warning(f"Failed with param {param}: {e}")
            # 診断情報を保存
            save_diagnostics(attempt, e)
    
    raise RuntimeError(f"All {max_retries} attempts failed: {last_error}")
```

### 新しい力場/水モデルの追加

`build_complex_system()`または`build_multi_ligand_system()`を拡張:

```python
# 力場マッピングに追加
forcefield_map = {
    "ff14SB": "leaprc.protein.ff14SB",
    "ff19SB": "leaprc.protein.ff19SB",  # 新規追加
    "ff14SBonlysc": "leaprc.protein.ff14SBonlysc",
}

# 水モデルマッピングに追加
water_source = {
    "tip3p": "leaprc.water.tip3p",
    "tip4pew": "leaprc.water.tip4pew",
    "opc": "leaprc.water.opc",
    "opc3": "leaprc.water.opc3",  # 新規追加
    "spce": "leaprc.water.spce",  # 新規追加
}

water_box = {
    "tip3p": "TIP3PBOX",
    "tip4pew": "TIP4PEWBOX",
    "opc": "OPCBOX",
    "opc3": "OPC3BOX",  # 新規追加
    "spce": "SPCBOX",   # 新規追加
}
```

### tleapスクリプトのカスタマイズ

複雑なシステム（膜タンパク質など）の場合、tleapスクリプトを直接編集:

```python
def build_membrane_system(protein_pdb, ligand_mol2, lipid_type="POPC"):
    """膜タンパク質システム構築の例"""
    
    leap_script = f"""
# 力場
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.lipid17  # 脂質力場

# リガンドパラメータ
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}

# タンパク質
protein = loadpdb {protein_pdb}

# 複合体作成
complex = combine {{protein LIG}}

# 膜の追加（PACKMOL等で事前準備が必要）
membrane = loadpdb membrane.pdb

# 結合
system = combine {{complex membrane}}

# 溶媒和（膜に合わせたボックス）
solvatebox system TIP3PBOX {{0 0 15}}

# イオン
addions system Na+ 0
addions system Cl- 0

saveamberparm system system.parm7 system.rst7
quit
"""
    # ... 以下実行 ...
```

### パス処理のベストプラクティス

```python
# ✓ 推奨: resolve()で絶対パス化
input_path = Path(input_file).resolve()
output_path = Path(output_dir).resolve() / "output.mol2"

# ✓ 推奨: ensure_directory()で作成
from common.utils import ensure_directory
ensure_directory(output_dir)

# ✗ 非推奨: 相対パスのまま使用
# input_path = Path(input_file)  # cwdに依存して失敗する可能性

# ✓ 推奨: 外部コマンドには絶対パスを渡す
wrapper.run(['-i', str(input_path.resolve()), '-o', str(output_path.resolve())])
```

### ログ出力のガイドライン

```python
from common.utils import setup_logger
logger = setup_logger(__name__)

# レベル別使い分け
logger.debug("詳細なデバッグ情報")
logger.info("処理開始: {file}")           # 通常の進捗
logger.warning("非致命的な問題: {msg}")   # 継続可能な警告
logger.error("エラー発生: {e}")           # 失敗時
```

### テストの書き方

`notebooks/test_amber_prep.ipynb`のパターンを参考に:

```python
# 関数取得ヘルパー
def get_callable(tool):
    """FunctionToolまたは通常の関数からcallableを取得"""
    if hasattr(tool, 'fn'):
        return tool.fn
    return tool

# テスト実行
my_tool = get_callable(amber_module.my_new_tool)
result = my_tool(input_file="test.pdb", option1="value")

# 検証
assert result['status'] == 'success'
assert Path(result['output_file']).exists()
```

### デバッグのヒント

1. **sqm失敗時**: `diagnostics/sqm_attempt*.out`を確認
2. **tleap失敗時**: `leap.log`を確認、特に"ERROR"行
3. **RDKit失敗時**: `Chem.MolFromPDBFile(..., sanitize=False)`で読み込んでから`Chem.SanitizeMol()`を別途実行
4. **パス問題**: `logger.info(f"Path: {path.resolve()}")`で絶対パス確認

---

## 拡張アイデア

### 1. 金属イオン対応
金属含有タンパク質のパラメータ化（MCPB.pyとの連携）

### 2. 共有結合リガンド
リガンドとタンパク質の共有結合処理

### 3. 自動プロトン化状態決定
PropKa/PDB2PQRとの統合

### 4. FEP/TI計算用セットアップ
自由エネルギー摂動計算用のラムダウィンドウ生成

### 5. コース粒子化
Martiniフォースフィールドへの変換

---

## 参考文献

- [Amber Manual](https://ambermd.org/Manuals.html)
- [OpenMM Documentation](https://openmm.org/documentation)
- [Boltz-2 GitHub](https://github.com/jwohlwend/boltz)
- [RDKit SMARTS Theory](https://www.rdkit.org/docs/RDKit_Book.html#smarts-support)
- [GAFF2 Force Field Paper](https://doi.org/10.1021/acs.jctc.5b00255)

