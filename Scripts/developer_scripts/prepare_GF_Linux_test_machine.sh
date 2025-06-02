#!/bin/bash
set -o errexit -o nounset -o pipefail

# Color output setup
if [ -t 1 ]; then
  COLOR_BOLD="\033[1m"
  COLOR_GREEN="\033[32m"
  COLOR_RESET="\033[0m"
else
  COLOR_BOLD=""
  COLOR_GREEN=""
  COLOR_RESET=""
fi

color_echo() {
  # Print with color unless line starts with two spaces
  if [[ "$1" =~ ^\ \  ]]; then
    echo "$1"
  else
    echo -e "${COLOR_BOLD}${COLOR_GREEN}$1${COLOR_RESET}"
  fi
}

color_printf() {
  # Print with color unless first arg starts with two spaces
  local fmt="$1"
  if [[ "$fmt" =~ ^\ \  ]]; then
    printf "$@"
  else
    echo -e "${COLOR_BOLD}${COLOR_GREEN}"
    printf "$@"
    echo -e "${COLOR_RESET}"
  fi
}

print_command() {
  if [ ! -t 0 ]; then
    while IFS= read -r line; do
      printf '  | %s\n' "$line" >&2
    done
  fi
  printf '  > %s\n' "$*" >&2
}

install_packages() {
  local packages=("$@")
  if [ "${#packages[@]}" -eq 0 ]; then
    return
  fi
  color_echo "Installing '${packages[*]}' (requires sudo privileges)..."
  local package_manager=""
  if command -v dnf &>/dev/null; then
    package_manager=dnf
  elif command -v apt &>/dev/null; then
    package_manager=apt
  elif command -v pacman &>/dev/null; then
    package_manager=pacman
  else
    echo "Could not find a supported package manager to install '${packages[*]}'."
    exit 1
  fi
  case "$package_manager" in
  dnf)
    $SUDO_OR_PRINT dnf install -y "${packages[@]}"
    ;;
  apt)
    $SUDO_OR_PRINT apt update
    $SUDO_OR_PRINT apt install -y "${packages[@]}"
    ;;
  pacman)
    $SUDO_OR_PRINT pacman -Sy --noconfirm "${packages[@]}"
    ;;
  esac
}

SUDO_OR_PRINT="sudo"
if [ "${1:-}" = "--dry-run" ]; then
  SUDO_OR_PRINT=print_command
  echo "Running in dry-run mode. No changes will be made."
fi

# Ensure gh CLI and jq are installed
required_pkgs=()
for pkg in gh jq; do
  command -v "$pkg" &>/dev/null || required_pkgs+=("$pkg")
done
install_packages "${required_pkgs[@]}"

if ! gh auth status &>/dev/null; then
  echo "GitHub CLI is not authenticated. Attempting to log in..."
  gh auth login || {
    echo "Failed to authenticate with GitHub CLI."
    exit 1
  }
fi

# Fetch members and real names using GraphQL
GF_GRAPHQL_QUERY='
query {
  organization(login: "CGAL") {
    team(slug: "GeometryFactory") {
      members(first: 100) {
        edges {
          node {
            login
            name
            publicKeys(first: 100) {
              edges {
                node {
                  id
                  key
                }
              }
            }
          }
        }
      }
    }
  }
}'

# Fetch team members and their names from GitHub GraphQL API
GRAPH_QL_RESULT=$(gh api graphql -f query="$GF_GRAPHQL_QUERY" --jq '.data.organization.team.members.edges.[].node')
if [ -z "$GRAPH_QL_RESULT" ]; then
  echo "Error: No members found in the GeometryFactory team." >&2
  exit 1
fi

GF_MEMBERS=()
declare -A GF_REAL_NAMES=()
declare -A GF_SSH_KEYS=()

# Parse the JSON result for logins, and names
while IFS=$'\t' read -r login name; do
  GF_MEMBERS+=("$login")
  GF_REAL_NAMES["$login"]="$name"
done < <(
  jq -r '[ .login, .name ] | @tsv' <<<"$GRAPH_QL_RESULT"
)

# Parse the SSH keys and associate them with the respective members
for login in "${GF_MEMBERS[@]}"; do
  keys=$(jq -r --arg login "$login" 'select(.login == $login) | (.publicKeys.edges[].node.key // empty)' <<<"$GRAPH_QL_RESULT")
  GF_SSH_KEYS["$login"]=""
  while IFS=$'\n' read -r key; do
    [[ -z "$key" ]] && continue
    GF_SSH_KEYS["$login"]+="$key"$'\n'
  done <<<"$keys"
done

# Mapping from GitHub login to local username (if different)
declare -A GF_USERNAMES=(
  ["janetournois"]="jtournois"
  ["LeoValque"]="lvalque"
  ["MaelRL"]="mrouxell"
  ["SaillantNicolas"]="nsaillant"
)

# Correcting the real names for specific users
GF_REAL_NAMES["LeoValque"]="Léo Valque"
GF_REAL_NAMES["MaelRL"]="Maël Rouxel-Labbé"
GF_REAL_NAMES["sloriot"]="Sébastien Loriot"

# Mapping from GitHub login to extra groups (if any)
declare -A GF_GROUPS=(
  ["lrineau"]="wheel,docker"
  ["SaillantNicolas"]="wheel,docker"
  ["sloriot"]="wheel,docker"
)

color_echo "Members of the geometryfactory team:"
for m in "${GF_MEMBERS[@]}"; do
  keys="${GF_SSH_KEYS[$m]}"
  nb_of_keys=$(printf '%s' "$keys" | wc -l)
  color_printf "  %s (%s) has %d ssh keys in Github\n" "$m" "${GF_REAL_NAMES[$m]:-$m}" "$nb_of_keys"
done

# Ensure groups exists
for group in geometryfactory docker; do
  if ! getent group "$group" >/dev/null; then
    $SUDO_OR_PRINT groupadd "$group"
  fi
done

for member in "${GF_MEMBERS[@]}"; do
  user="${GF_USERNAMES[$member]:-$member}"
  groups="geometryfactory"
  if [[ -n "${GF_GROUPS[$member]:-}" ]]; then
    groups="geometryfactory,${GF_GROUPS[$member]}"
  fi

  if ! getent passwd "$user" >/dev/null; then
    color_echo "Creating user $user and adding to group(s) $groups."
    $SUDO_OR_PRINT useradd -m -s /bin/bash -c "${GF_REAL_NAMES[$member]:-$user}" "$user" -G "$groups"
  else
    color_echo "User $user already exists, adding to group(s) $groups."
    $SUDO_OR_PRINT usermod -aG "$groups" -c "${GF_REAL_NAMES[$member]:-$user}" "$user"
  fi

  $SUDO_OR_PRINT chmod o+rx "/home/$user"

  $SUDO_OR_PRINT mkdir -p "/home/$user/.ssh"
  $SUDO_OR_PRINT chown "$user:$user" "/home/$user/.ssh"
  $SUDO_OR_PRINT chmod 700 "/home/$user/.ssh"

  # Add SSH keys from GraphQL
  color_echo "  Adding SSH keys for user $member from Github."
  count=0
  auth_keys_file="/home/$user/.ssh/authorized_keys"
  if ! sudo test -f "$auth_keys_file"; then
    $SUDO_OR_PRINT install -m 600 -o "$user" -g "$user" /dev/null "$auth_keys_file"
  fi
  if [[ -n "${GF_SSH_KEYS[$member]:-}" ]]; then
    while IFS=$'\n' read -r key; do
      [[ -z "$key" ]] && continue
      if ! $SUDO_OR_PRINT grep -qF -- "$key" "$auth_keys_file" </dev/null; then
        if [ "$SUDO_OR_PRINT" = "sudo" ]; then
          echo "$key" | sudo tee -a "$auth_keys_file" >/dev/null
        fi
        count=$((count + 1))
      fi
    done <<<"${GF_SSH_KEYS[$member]}"
  fi
  color_echo "  Added $count SSH keys for user $user."
  key_count=$(sudo bash -c "[ -f '$auth_keys_file' ] && wc -l <'$auth_keys_file' || echo 0")
  color_echo "  User $user has $key_count SSH keys in $auth_keys_file"
done

# Add cgaltest user if not exists, and create SSH key
if ! getent passwd cgaltest >/dev/null; then
  color_echo "Creating user cgaltest."
  $SUDO_OR_PRINT useradd -m -s /bin/bash -c "CGAL Test User" cgaltest
fi

$SUDO_OR_PRINT mkdir -p /home/cgaltest/.ssh
$SUDO_OR_PRINT chown cgaltest:cgaltest /home/cgaltest/.ssh
$SUDO_OR_PRINT chmod 700 /home/cgaltest/.ssh
$SUDO_OR_PRINT usermod -aG docker cgaltest

keyfile="/home/cgaltest/.ssh/id_ed25519"
if ! $SUDO_OR_PRINT test -f "$keyfile"; then
  color_echo "Generating SSH key for cgaltest."
  if [ "$SUDO_OR_PRINT" = "sudo" ]; then
    sudo -u cgaltest ssh-keygen -t ed25519 -N "" -f "$keyfile"
  else
    ssh-keygen -t ed25519 -N "" -f "$keyfile"
  fi
else
  color_echo "SSH key for cgaltest already exists."
  $SUDO_OR_PRINT cat "$keyfile.pub"
fi

color_echo "All users have been processed."

color_echo "Remove users not in the GeometryFactory team."
current_members=$(getent group geometryfactory | awk -F: '{print $4}' | tr ',' ' ')
for user in $current_members; do
  # Map local usernames to GitHub logins if needed
  login="$user"
  for k in "${!GF_USERNAMES[@]}"; do
    if [[ "${GF_USERNAMES[$k]}" == "$user" ]]; then
      login="$k"
      break
    fi
  done
  found=false
  for m in "${GF_MEMBERS[@]}"; do
    if [[ "$login" == "$m" ]]; then
      found=true
      break
    fi
  done
  if ! $found; then
    color_echo "  Removing user $user but keep its home directory."
    read -p "  Are you sure you want to remove user $user? [y/N] " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
      color_echo "  Skipping removal of user $user."
      continue
    fi
    $SUDO_OR_PRINT userdel -f "$user"
  fi
done

# Optionally install bat, strace, and podman-docker if not present
need_install=()
for cmd in bat strace; do
  command -v "$cmd" &>/dev/null || need_install+=("$cmd")
done

# Check if docker is missing or not emulated by podman
need_podman_docker=false
if ! command -v docker &>/dev/null; then
  need_podman_docker=true
elif docker --version 2>&1 | grep -q "podman"; then
  need_podman_docker=true
fi

if $need_podman_docker; then
  if [[ ! " ${need_install[*]} " =~ " podman-docker " ]]; then
    need_install+=("podman-docker")
  fi
fi

install_packages "${need_install[@]}"

color_echo "Update the systemd tmpfiles configuration for podman."

$SUDO_OR_PRINT touch /etc/containers/nodocker

$SUDO_OR_PRINT mkdir -p /etc/tmpfiles.d
sed 's|podman 0700 root root|podman 0700 root docker|g' /usr/lib/tmpfiles.d/podman.conf |
  $SUDO_OR_PRINT tee /etc/tmpfiles.d/podman.conf >/dev/null

color_echo "Update the systemd socket configuration for podman."

$SUDO_OR_PRINT mkdir -p /etc/systemd/system/podman.socket.d
$SUDO_OR_PRINT tee /etc/systemd/system/podman.socket.d/override.conf >/dev/null <<EOF
[Socket]
SocketMode=0660
SocketUser=root
SocketGroup=docker
SELinuxContextFromNet=true
SELinuxContext=system_u:object_r:container_var_run_t:s0
EOF

color_echo "Reloading systemd and enabling podman.socket."
$SUDO_OR_PRINT systemctl daemon-reload
$SUDO_OR_PRINT systemctl enable podman.socket
$SUDO_OR_PRINT systemctl restart podman.socket

color_echo "All done!"
