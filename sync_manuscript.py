# -*- coding: utf-8 -*-
"""
一键同步脚本：你手动改完 tsc-gnn-repo/gate1/rewiring_study/manuscript_v6_humanized.md 后，
双击同目录的 sync_manually.bat 即可：

  1) 把你在【仓库内】改好的 md 同步到 build 上下文（父目录那份，build_docx.py 读取它）；
  2) 用受控 venv（python-docx）重新生成 docx；
  3) 把新 docx 写回仓库 + Patterns 投稿包；
  4) git 提交并推送（HTTP/1.1 绕过公司代理 reset）。

约定：仓库里的 md 是「编辑源」，父目录那份只是 build 镜像（会被覆盖）。
所有路径都相对本脚本位置推导，不依赖工作目录，也不硬编码中文路径。
"""
import os
import sys
import shutil
import datetime
import subprocess

# ---------- 路径推导（全部相对脚本自身） ----------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))   # = tsc-gnn-repo/
REPO_ROOT  = SCRIPT_DIR
PARENT     = os.path.dirname(REPO_ROOT)                   # = 科研/虚拟敲除/
REW        = os.path.join(REPO_ROOT, "gate1", "rewiring_study")
PARENT_REW = os.path.join(PARENT, "gate1", "rewiring_study")

REPO_MD      = os.path.join(REW, "manuscript_v6_humanized.md")
PARENT_MD    = os.path.join(PARENT_REW, "manuscript_v6_humanized.md")
PARENT_DOCX  = os.path.join(PARENT_REW, "manuscript_v6_humanized.docx")
REPO_DOCX    = os.path.join(REW, "manuscript_v6_humanized.docx")
SUB_PKG_DOCX = os.path.join(PARENT_REW, "Patterns_submission_package", "Manuscript_TSC-GNN_v6.docx")
BUILD_SCRIPT = os.path.join(PARENT, "build_docx.py")

# venv（受控，python-docx 1.2.0 已装）；缺失时回退到 managed python 3.13.12 自建
MANAGED_PY = r"C:/Users/lm962/.workbuddy/binaries/python/versions/3.13.12/python.exe"
VENV_PY    = r"C:/Users/lm962/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
GH         = r"C:/Users/lm962/AppData/Local/gh-cli/bin/gh.exe"


def run(cmd, **kw):
    print(">>", " ".join(cmd))
    return subprocess.run(cmd, **kw)


def ensure_venv():
    """确保存在可用的 venv 且装了 python-docx。"""
    if not os.path.exists(VENV_PY):
        print("[setup] 创建 venv ...")
        subprocess.run([MANAGED_PY, "-m", "venv",
                        os.path.dirname(os.path.dirname(VENV_PY))], check=True)
    try:
        subprocess.run([VENV_PY, "-c", "import docx"],
                       check=True, capture_output=True)
    except subprocess.CalledProcessError:
        print("[setup] 安装 python-docx ...")
        subprocess.run([VENV_PY, "-m", "pip", "install", "python-docx"], check=True)
    return VENV_PY


def main():
    # 0) 仓库内 md 必须存在
    if not os.path.exists(REPO_MD):
        print("ERROR: 找不到编辑源", REPO_MD)
        sys.exit(1)

    # 1) 同步 md → build 上下文（父目录那份）
    shutil.copy2(REPO_MD, PARENT_MD)
    print("[1] 已同步 md → build 上下文")

    # 2) 重建 docx
    built = False
    if os.path.exists(BUILD_SCRIPT):
        try:
            ensure_venv()
            r = run([VENV_PY, BUILD_SCRIPT])
            built = (r.returncode == 0)
        except Exception as e:
            print("[warn] docx 重建异常：", e)
        if not built:
            print("[warn] docx 重建失败，本次仅同步 md")
    else:
        print("[warn] 找不到 build_docx.py，跳过 docx 重建")

    # 3) 写回 docx（仓库 + 投稿包）
    if built and os.path.exists(PARENT_DOCX):
        shutil.copy2(PARENT_DOCX, REPO_DOCX)
        os.makedirs(os.path.dirname(SUB_PKG_DOCX), exist_ok=True)
        shutil.copy2(PARENT_DOCX, SUB_PKG_DOCX)
        print("[2] 已重建并写回 docx（仓库 + 投稿包）")
    elif os.path.exists(REPO_DOCX):
        # build 失败但仓库里已有旧 docx：至少把 md 改动带进去，docx 维持现状
        print("[2] docx 未重建，保持仓库现有 docx")

    # 4) git 提交 + 推送
    os.chdir(REPO_ROOT)
    st = subprocess.run(["git", "status", "--porcelain"],
                        capture_output=True, text=True)
    if not st.stdout.strip():
        print("[3] 无改动，无需同步。")
        return

    today = datetime.date.today().isoformat()
    # 只提交与稿件同步相关的文件（避免误带其它改动）
    add_targets = [
        "gate1/rewiring_study/manuscript_v6_humanized.md",
        "gate1/rewiring_study/manuscript_v6_humanized.docx",
        "sync_manuscript.py",
        "sync_manually.bat",
    ]
    add_list = [t for t in add_targets if os.path.exists(os.path.join(REPO_ROOT, t))]
    if not add_list:
        print("[3] 无可提交的目标文件，跳过。")
        return
    run(["git", "add"] + add_list)
    run(["git", "commit", "-m", f"manuscript update {today} (md+docx synced)"])

    # 推送前检查 gh 登录态
    authed = os.path.exists(GH) and (
        subprocess.run([GH, "auth", "status"],
                       capture_output=True, text=True).returncode == 0
    )
    if not authed:
        print("[push] gh 尚未登录，已提交到本地。请在本机终端执行：")
        print(f'       "{GH}" auth login')
        print(f'       "{GH}" auth setup-git')
        print("       然后重跑本脚本，或手动推送：")
        print("       git -c http.version=HTTP/1.1 -c http.postBuffer=524288000 push origin main")
        return

    run(["git", "-c", "http.version=HTTP/1.1",
         "-c", "http.postBuffer=524288000", "push", "origin", "main"])
    print("[done] 已提交并推送。")


if __name__ == "__main__":
    main()
