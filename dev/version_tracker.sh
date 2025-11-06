#!/usr/bin/env bash

# version_tracker.sh
# Description: increments version, records changes per user and globally, tracks activity in CSV.

set -euo pipefail
IFS=$'\n\t'

# -------------------- #
# Helper: increment_version
# -------------------- #
increment_version() {
    local version="$1"
    local part="${2:-z}"  # default to patch (z)
    IFS='.' read -r major minor patch <<< "$version"

    case "$part" in
        x) major=$((major + 1)); minor=0; patch=0 ;;
        y) minor=$((minor + 1)); patch=0 ;;
        z) patch=$((patch + 1)) ;;
        *) echo "❌ Invalid part: $part. Use 'x', 'y', or 'z'." >&2; exit 1 ;;
    esac
    echo "$major.$minor.$patch"
}

# -------------------- #
# Setup
# -------------------- #
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

VERSION_FILE="${LOG_DIR}/VERSION"
GLOBAL_CHANGELOG="${LOG_DIR}/changelog.md"
USER_NAME="${USER:-$(whoami)}"
USER_LOG="${LOG_DIR}/${USER_NAME}.log"
ACTIVITY_CSV="${LOG_DIR}/activity.csv"

# -------------------- #
# Initialize files if missing
# -------------------- #
if [[ ! -f "$VERSION_FILE" ]]; then
    echo "0.0.0" > "$VERSION_FILE"
fi

if [[ ! -f "$ACTIVITY_CSV" ]]; then
    echo "timestamp,user,new_version" > "$ACTIVITY_CSV"
fi

current_version=$(<"$VERSION_FILE")

# -------------------- #
# Ask user which part to increment
# -------------------- #
read -p "-- Current version is $current_version. Increment which part? (x/y/z, default=z): " part
part="${part:-z}"

if [[ "$part" != "x" && "$part" != "y" && "$part" != "z" ]]; then
    echo "❌ Invalid input: $part. Please use 'x', 'y', or 'z'."
    exit 1
fi

new_version=$(increment_version "$current_version" "$part")

# -------------------- #
# Ask user for changelog entries
# -------------------- #
echo "-- Describe what you've done for version $new_version."
echo "/!\\ Enter each change on a separate line."
echo "/!\\ Press Ctrl+C and re-run the script to start over "
echo "/!\\ Press Enter on an empty line when done:"
entries=()
while IFS= read -r line; do
    [[ -z "$line" ]] && break
    entries+=("$line")
done

if [[ ${#entries[@]} -eq 0 ]]; then
    echo "❌ No changes provided. Exiting."
    exit 1
fi

# -------------------- #
# Update VERSION file
# -------------------- #
echo "$new_version" > "$VERSION_FILE"
echo "✅ Updated VERSION → $new_version"

# -------------------- #
# Write global changelog
# -------------------- #
{
    echo "## [$new_version] - $(date +'%Y-%m-%d') by @$USER_NAME"
    for e in "${entries[@]}"; do
        echo "- $e"
    done
    echo
} >> "$GLOBAL_CHANGELOG"

# -------------------- #
# Write per-user log
# -------------------- #
{
    echo "### [$new_version] - $(date +'%Y-%m-%d %H:%M:%S')"
    for e in "${entries[@]}"; do
        echo "- $e"
    done
    echo
} >> "$USER_LOG"

# -------------------- #
# Update activity CSV
# -------------------- #
timestamp=$(date +'%Y-%m-%d %H:%M:%S')
echo "${timestamp},${USER_NAME},${new_version}" >> "$ACTIVITY_CSV"

echo "📘 Logged in:"
echo "  - Global changelog: $GLOBAL_CHANGELOG"
echo "  - User log: $USER_LOG"
echo "  - Activity CSV: $ACTIVITY_CSV"
echo "✨ Version updated successfully from $current_version → $new_version."
